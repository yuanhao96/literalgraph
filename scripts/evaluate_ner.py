from seqeval.metrics import (
    classification_report,
    precision_score,
    recall_score,
    f1_score
)
from seqeval.metrics.sequence_labeling import get_entities
from tqdm import tqdm
import re
import pandas as pd
from datasets import load_dataset
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from information_extraction.NER import BiomedicalNER

# Your valid model entity types
model_entity_types = {
    'gene', 'disease', 'chemical', 'organism',
    'anatomical', 'cell_type', 'cell_line', 'variant'
}

# Mapping from GENIA fine-grained types to model types
genia_to_model_map = {
    'protein': 'gene', 
    # 'protein_N/A': 'gene', 
    'protein_molecule': 'gene',
    # 'protein_subunit': 'gene', 
    'protein_complex': 'gene', 
    'protein_family_or_group': 'gene',
    # 'protein_domain_or_region': 'gene', 
    # 'protein_substructure': 'gene',
    # 'DNA_molecule': 'gene',
    # 'DNA_domain_or_region': 'gene', 
    'DNA_family_or_group': 'gene',
    # 'DNA_substructure': 'gene', 
    # 'polynucleotide': 'gene', 'nucleotide': 'gene',
    # 'RNA_molecule': 'gene', 'RNA_family_or_group': 'gene', 'RNA_substructure': 'gene',
    'cell_type': 'cell_type', 'mono_cell': 'cell_type', 'multi_cell': 'cell_type',
    'cell_line': 'cell_line',
    # 'body_part': 'anatomical', 'tissue': 'anatomical', 'cell_component': 'anatomical',
    'virus': 'organism', 'organism': 'organism',
    # 'peptide': 'chemical', 
    'lipid': 'chemical', 
    'amino_acid_monomer': 'chemical',
    'other_organic_compound': 'chemical', 'inorganic': 'chemical', 'atom': 'chemical',
    # 'variant': 'variant',
    # 'other_artificial_source': 'chemical',
    # 'other_name': 'gene',
    # '(NEITHER_NOR RNA_molecule RNA_molecule)': 'gene'
}

def convert_genia_label(label):
    if label == "O":
        return "O"
    match = re.match(r'^[BI]-(.+)', label)
    if not match:
        return "O"
    base_type = match.group(1)
    mapped_type = genia_to_model_map.get(base_type)
    if mapped_type in model_entity_types:
        return label[0] + '-' + mapped_type
    else:
        return "O"

def evaluate_biomedical_ner_on_genia(ner, dataset, evaluate_with_llm=False):
    true_labels = []
    pred_labels = []

    for example in tqdm(dataset, desc="Evaluating"):
        tokens = example["tokens"]
        raw_labels = example["labels"]
        text = " ".join(tokens)

        # Convert GENIA fine-grained labels to model-compatible tags
        gold_labels = [convert_genia_label(lbl) for lbl in raw_labels]
        true_labels.append(gold_labels)

        # Predict using your model
        result = ner.extract_entities(text, entity_types=['gene', 'chemical', 'organism', 'anatomical', 'cell_type', 'cell_line'], confidence_threshold=0.3, ground_entities=False, evaluate_with_llm=evaluate_with_llm)
        predicted_tags = ["O"] * len(tokens)

        # Char → token mapping
        char_to_token = {}
        cur_char = 0
        for idx, token in enumerate(tokens):
            for c in range(len(token)):
                char_to_token[cur_char + c] = idx
            cur_char += len(token) + 1  # +1 for space

        for entity_type, entities in result.items():
            if entity_type not in model_entity_types:
                continue
            for ent in entities:
                start_token = char_to_token.get(ent['start'])
                end_token = char_to_token.get(ent['end'] - 1)
                if start_token is None or end_token is None:
                    continue
                for i in range(start_token, end_token + 1):
                    tag = f"{'B' if i == start_token else 'I'}-{entity_type}"
                    if i < len(predicted_tags):
                        predicted_tags[i] = tag

        pred_labels.append(predicted_tags)

    # Compute metrics
    global_precision = precision_score(true_labels, pred_labels)
    global_recall = recall_score(true_labels, pred_labels)
    global_f1 = f1_score(true_labels, pred_labels)

    print("\n🔹 Global (micro-averaged) scores:")
    print(f"Precision: {global_precision:.4f}")
    print(f"Recall:    {global_recall:.4f}")
    print(f"F1-score:  {global_f1:.4f}")

    print("\n🔹 Per-entity-type classification report:")
    label_report = classification_report(true_labels, pred_labels, output_dict=True, digits=4)
    df = pd.DataFrame(label_report).transpose()
    display_cols = ['precision', 'recall', 'f1-score', 'support']
    df = df.loc[[k for k in df.index if k not in ['micro avg', 'macro avg', 'weighted avg', 'accuracy']], display_cols]
    print(df.to_string())
    
    return {
        'global': {
            'precision': global_precision,
            'recall': global_recall,
            'f1': global_f1
        },
        'per_label': df
    }

def evaluate_biomedical_ner_on_ncbi_disease(ner, dataset, evaluate_with_llm=False):
    # Map tag indices to names (e.g., 0->'O', 1->'B-Disease', etc.)
    label_list = dataset.features['ner_tags'].feature.names
    
    true_labels = []
    pred_labels = []

    for example in tqdm(dataset, desc="Evaluating NCBI Disease"):
        tokens = example['tokens']
        gold_ids = example['ner_tags']
        gold_tags = [label_list[i] for i in gold_ids]
        true_labels.append(gold_tags)

        text = " ".join(tokens)

        # Predict using the 'disease' model
        result = ner.extract_entities(
            text, entity_types=['disease'], confidence_threshold=0.3, ground_entities=False, evaluate_with_llm=evaluate_with_llm
        )
        predicted_tags = ["O"] * len(tokens)

        # Map character positions to token indices
        char_to_token = {}
        cur_char = 0
        for idx, tok in enumerate(tokens):
            for c in range(len(tok)):
                char_to_token[cur_char + c] = idx
            cur_char += len(tok) + 1

        for ent in result.get('disease', []):
            start_token = char_to_token.get(ent['start'])
            end_token = char_to_token.get(ent['end'] - 1)
            if start_token is None or end_token is None:
                continue
            for i in range(start_token, end_token + 1):
                tag = "B-Disease" if i == start_token else "I-Disease"
                if i < len(predicted_tags):
                    predicted_tags[i] = tag

        pred_labels.append(predicted_tags)

    # Compute metrics
    global_precision = precision_score(true_labels, pred_labels)
    global_recall = recall_score(true_labels, pred_labels)
    global_f1 = f1_score(true_labels, pred_labels)

    print("\n🔹 Global (micro-averaged) scores:")
    print(f"Precision: {global_precision:.4f}")
    print(f"Recall:    {global_recall:.4f}")
    print(f"F1-score:  {global_f1:.4f}")

    print("\n🔹 Per-label classification report:")
    print(classification_report(true_labels, pred_labels, digits=4))

    return {
        'global': {
            'precision': global_precision,
            'recall': global_recall,
            'f1': global_f1
        }
    }

def evaluate_biomedical_ner_on_bc2gm(ner, dataset, evaluate_with_llm=False):
    # Map tag indices to names (e.g., 0->'O', 1->'B-Disease', etc.)
    label_list = dataset.features['ner_tags'].feature.names
    
    true_labels = []
    pred_labels = []

    for example in tqdm(dataset, desc="Evaluating BC2GM"):
        tokens = example['tokens']
        gold_ids = example['ner_tags']
        gold_tags = [label_list[i] for i in gold_ids]
        true_labels.append(gold_tags)

        text = " ".join(tokens)
        # print(text)

        # Predict using the 'gene' model
        result = ner.extract_entities(
            text, entity_types=['gene'], confidence_threshold=0.3, ground_entities=False, evaluate_with_llm=evaluate_with_llm
        )
        predicted_tags = ["O"] * len(tokens)
        # print(result)

        # Map character positions to token indices
        char_to_token = {}
        cur_char = 0
        for idx, tok in enumerate(tokens):
            for c in range(len(tok)):
                char_to_token[cur_char + c] = idx
            cur_char += len(tok) + 1

        for ent in result.get('gene', []):
            start_token = char_to_token.get(ent['start'])
            end_token = char_to_token.get(ent['end'] - 1)
            if start_token is None or end_token is None:
                continue
            for i in range(start_token, end_token + 1):
                tag = "B-GENE" if i == start_token else "I-GENE"
                if i < len(predicted_tags):
                    predicted_tags[i] = tag

        pred_labels.append(predicted_tags)
        # print(gold_tags)
        # print(predicted_tags)

    # Compute metrics
    global_precision = precision_score(true_labels, pred_labels)
    global_recall = recall_score(true_labels, pred_labels)
    global_f1 = f1_score(true_labels, pred_labels)

    print("\n🔹 Global (micro-averaged) scores:")
    print(f"Precision: {global_precision:.4f}")
    print(f"Recall:    {global_recall:.4f}")
    print(f"F1-score:  {global_f1:.4f}")

    print("\n🔹 Per-label classification report:")
    print(classification_report(true_labels, pred_labels, digits=4))

    return {
        'global': {
            'precision': global_precision,
            'recall': global_recall,
            'f1': global_f1
        }
    }

bionlp_to_model_map = {
    'PROTEIN': 'gene',
    'ORGANISM': 'organism',
    'CHEMICAL': 'chemical'
}

def convert_bionlp_label(label):
    if label == "O":
        return "O"
    match = re.match(r'^([BI])-(.+)$', label)
    if not match:
        return "O"
    prefix, fine_type = match.groups()
    coarse = bionlp_to_model_map.get(fine_type)
    if coarse in model_entity_types:
        return f"{prefix}-{coarse}"
    return "O"

def evaluate_biomedical_ner_on_bionlp11id(ner, dataset, evaluate_with_llm=False):
    true_labels = []
    pred_labels = []

    for example in tqdm(dataset, desc="Evaluating BioNLP11ID-ggp"):
        tokens = example['tokens']
        raw_tags = example['ner_tags']
        gold_tags = [convert_bionlp_label(tag) for tag in raw_tags]

        text = " ".join(tokens)
        predicted_tags = ["O"] * len(tokens)

        # Char → token index map
        char_to_token = {}
        cur_char = 0
        for idx, tok in enumerate(tokens):
            for c in range(len(tok)):
                char_to_token[cur_char + c] = idx
            cur_char += len(tok) + 1

        result = ner.extract_entities(text, entity_types=['gene', 'chemical', 'organism'], confidence_threshold=0.0, ground_entities=False, evaluate_with_llm=evaluate_with_llm)

        for entity_type, entities in result.items():
            if entity_type not in model_entity_types:
                continue
            for ent in entities:
                start_token = char_to_token.get(ent['start'])
                end_token = char_to_token.get(ent['end'] - 1)
                if start_token is None or end_token is None:
                    continue
                for i in range(start_token, end_token + 1):
                    tag = f"{'B' if i == start_token else 'I'}-{entity_type}"
                    if i < len(predicted_tags):
                        predicted_tags[i] = tag

        if len(predicted_tags) != len(gold_tags):
            print(f"[Warning] Length mismatch for id {example.get('id', 'unknown')} — skipping")
            continue

        true_labels.append(gold_tags)
        pred_labels.append(predicted_tags)

    # Metrics
    global_precision = precision_score(true_labels, pred_labels)
    global_recall = recall_score(true_labels, pred_labels)
    global_f1 = f1_score(true_labels, pred_labels)

    print("\n🔹 Global (micro-averaged) scores:")
    print(f"Precision: {global_precision:.4f}")
    print(f"Recall:    {global_recall:.4f}")
    print(f"F1-score:  {global_f1:.4f}")

    print("\n🔹 Per-label classification report:")
    report = classification_report(true_labels, pred_labels, output_dict=True, digits=4)
    df = pd.DataFrame(report).transpose()
    df = df.loc[[k for k in df.index if k not in ['micro avg', 'macro avg', 'weighted avg', 'accuracy']], ['precision', 'recall', 'f1-score', 'support']]
    print(df.to_string())

    return {
        'global': {
            'precision': global_precision,
            'recall': global_recall,
            'f1': global_f1
        },
        'per_label': df
    }

if __name__ == "__main__":
    ner = BiomedicalNER()
    # genia_test = load_dataset("enoriega/GENIA-Term-Corpus", split="test")
    # genia_test = genia_test.select(range(1000))
    # disease_test = load_dataset("ncbi/ncbi_disease", split="test")
    # gene_test = load_dataset("omniquad/BioNLP11ID-ggp-IOB", split="test")
    # gene_test = gene_test.select(range(1000))
    gene_test = load_dataset("omniquad/BioNLP11ID-ggp-IOB", split="test")
    gene_test = gene_test.select(range(1000))

    # metrics = evaluate_biomedical_ner_on_genia(ner, genia_test, evaluate_with_llm=True)
    # metrics = evaluate_biomedical_ner_on_ncbi_disease(ner, disease_test, evaluate_with_llm=True)
    # metrics = evaluate_biomedical_ner_on_bc2gm(ner, gene_test, evaluate_with_llm=True)
    metrics = evaluate_biomedical_ner_on_bionlp11id(ner, gene_test, evaluate_with_llm=True)