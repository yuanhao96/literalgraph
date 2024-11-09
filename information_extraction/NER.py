from transformers import AutoTokenizer, AutoModelForTokenClassification
from transformers import pipeline
import torch
from entity_mapping.gilda_grounders import (
    Gene_Grounder, Disease_Grounder, Chemical_Grounder, 
    Organism_Grounder, Anatomy_Grounder, Variant_Grounder
)

class BiomedicalNER:
    def __init__(self, use_local_grounders=True):
        # Initialize tokenizers and models for all entity types
        self.models = {
            'gene': "pruas/BENT-PubMedBERT-NER-Gene",
            'disease': "pruas/BENT-PubMedBERT-NER-Disease",
            'cell_type': "pruas/BENT-PubMedBERT-NER-Cell-Type",
            'anatomical': "pruas/BENT-PubMedBERT-NER-Anatomical",
            'variant': "pruas/BENT-PubMedBERT-NER-Variant",
            'cell_line': "pruas/BENT-PubMedBERT-NER-Cell-Line",
            'chemical': "pruas/BENT-PubMedBERT-NER-Chemical",
            'organism': "pruas/BENT-PubMedBERT-NER-Organism"
        }
        
        # Initialize grounders
        if use_local_grounders:
            self.grounders = {
                'gene': Gene_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/gene.json'),
                'disease': Disease_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/disease.json'),
                'chemical': Chemical_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/chemical.json'),
                'organism': Organism_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/organism.json'),
                'anatomical': Anatomy_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/anatomy.json'),
                'cell_type': Anatomy_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/anatomy.json'),
                'cell_line': Anatomy_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/anatomy.json'),
                'variant': Variant_Grounder(),
            }
        else:
            self.grounders = {
                'gene': Gene_Grounder(),
                'disease': Disease_Grounder(),
                'chemical': Chemical_Grounder(),
                'organism': Organism_Grounder(),
                'anatomical': Anatomy_Grounder(),
                'cell_type': Anatomy_Grounder(),
                'cell_line': Anatomy_Grounder(),
                'variant': Variant_Grounder(),
            }
        
        # Initialize NER pipelines
        self.pipelines = {}
        for entity_type, model_name in self.models.items():
            try:
                tokenizer = AutoTokenizer.from_pretrained(model_name)
                model = AutoModelForTokenClassification.from_pretrained(model_name)
                self.pipelines[entity_type] = pipeline("ner", model=model, tokenizer=tokenizer)
            except Exception as e:
                print(f"Warning: Failed to load {entity_type} model: {str(e)}")

    def ground_entity(self, text, entity_type):
        """
        Ground an entity using the appropriate grounder
        """
        if entity_type not in self.grounders:
            return []
        
        terms = self.grounders[entity_type].ground(text)
        return [(m.term.get_curie(), m.score) for m in terms]

    def extract_entities(self, text, entity_types='all', confidence_threshold=0.5, ground_entities=True):
        """
        Extract and optionally ground entities from input text
        Args:
            text (str): Input text to analyze
            entity_types (str or list): Type(s) of entities to extract
            confidence_threshold (float): Minimum confidence score threshold
            ground_entities (bool): Whether to perform entity grounding
        Returns:
            dict: Dictionary containing entity mentions with their positions and groundings
        """
        results = {}
        
        # Determine which entity types to extract
        if entity_types == 'all':
            types_to_extract = self.pipelines.keys()
        elif isinstance(entity_types, str):
            types_to_extract = [entity_types]
        else:
            types_to_extract = entity_types
            
        for entity_type in types_to_extract:
            if entity_type not in self.pipelines:
                print(f"Warning: Entity type '{entity_type}' not available")
                continue
                
            try:
                ner_results = self.pipelines[entity_type](text)
                entities = []
                
                filtered_results = [r for r in ner_results if r['score'] > confidence_threshold]
                
                i = 0
                while i < len(filtered_results):
                    current = filtered_results[i]
                    entity = current['word']
                    start = current['start']
                    end = current['end']
                    scores = [current['score']]
                    
                    j = i + 1
                    while j < len(filtered_results):
                        next_token = filtered_results[j]
                        if next_token['start'] == end or next_token['start'] == end + 1:
                            entity = text[start:next_token['end']]
                            end = next_token['end']
                            scores.append(next_token['score'])
                            j += 1
                        else:
                            break
                    
                    entity_info = {
                        'entity': entity,
                        'start': start,
                        'end': end,
                        'score': round(sum(scores) / len(scores), 3)
                    }
                    
                    # Add grounding information if requested
                    if ground_entities:
                        groundings = self.ground_entity(entity, entity_type)
                        if groundings:
                            entity_info['groundings'] = groundings
                    
                    entities.append(entity_info)
                    i = j if j > i + 1 else i + 1
                    
                results[entity_type] = entities
                
            except Exception as e:
                print(f"Error extracting {entity_type} entities: {str(e)}")
                
        return results

# Example usage
if __name__ == "__main__":
    ner = BiomedicalNER()
    
    text = """BRCA1 mutations and rs12345 were studied in MCF-7 breast cancer cells.
              The variant rs987654 is associated with drug response.
              The EGFR T790M variant was treated with gefitinib."""
    
    results = ner.extract_entities(text, entity_types='all', ground_entities=True)
    
    for entity_type, entities in results.items():
        if entities:
            print(f"\nDetected {entity_type.replace('_', ' ').title()}:")
            for entity in entities:
                print(f"Entity: {entity['entity']}")
                print(f"Position: {entity['start']}-{entity['end']}")
                print(f"Confidence: {entity['score']}")
                if 'groundings' in entity:
                    print("Groundings:")
                    for curie, score in entity['groundings']:
                        print(f"  - {curie} (score: {score:.3f})")
                print()
