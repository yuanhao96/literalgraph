# MedlineKG Source Code

This directory contains the source code for MedlineKG, a biomedical knowledge graph construction and evaluation system that processes PubMed literature and integrates multiple biomedical data sources.

## Project Overview

MedlineKG is designed to build a comprehensive biomedical knowledge graph by:
- Extracting entities from PubMed articles using state-of-the-art NER models
- Integrating data from multiple biomedical databases (dbSNP, Reactome, GO, PrimeKG, GWAS, etc.)
- Evaluating entity extraction performance using LLM-based validation
- Summarizing relationships between biomedical entities

## Directory Structure

### `/adapters/`
Contains BioCypher adapters for integrating various biomedical data sources:

- **`pubmed_adapter.py`** - Processes PubMed XML files to extract articles, sentences, and metadata
- **`journal_adapter.py`** - Handles journal information and metadata
- **`vocab_adapter.py`** - Manages vocabulary and ontology mappings
- **`dbsnp_adapter.py`** - Integrates dbSNP variant data
- **`reactome_adapter.py`** - Processes Reactome pathway data
- **`go_adapter.py`** - Integrates Gene Ontology annotations
- **`primekg_adapter.py`** - Incorporates PrimeKG knowledge graph data
- **`gwas_adapter.py`** - Processes GWAS association data
- **`ner_adapter.py`** - Handles named entity recognition results
- **`nodes.py`** - Defines node types and structures
- **`edges.py`** - Defines edge types and relationships

### `/information_extraction/`
Core NER and information extraction modules:

- **`NER.py`** - Main biomedical NER system using BERT-based models
  - Supports 8 entity types: gene, disease, chemical, organism, anatomical, cell_type, cell_line, variant
  - Uses PubMedBERT-based models from Hugging Face
  - Includes entity grounding with GILDA
  - LLM-based validation for low-confidence extractions
- **`gilda_grounders.py`** - Entity grounding utilities using GILDA
- **`relation_summarization.py`** - LLM-based relationship summarization between entities
- **`env.sh`** - Environment configuration script

### `/llm_evaluation/`
LLM-based evaluation and validation tools:

- **`llm_evaluator.py`** - Evaluates NER performance using LLM validation
  - Connects to Neo4j database
  - Validates entity extractions against ground truth
  - Calculates false positive rates

### `/scripts/`
Main execution scripts:

- **`build_kg.py`** - Main script to build the knowledge graph
  - Orchestrates all adapters to process data sources
  - Uses BioCypher for graph construction
  - Handles batch processing and error management
- **`evaluate_ner.py`** - Comprehensive NER evaluation script
  - Evaluates on multiple datasets (GENIA, NCBI Disease, BC2GM, BioNLP11ID)
  - Supports both traditional and LLM-based evaluation
  - Generates detailed performance metrics

### `/utils/`
Utility functions and helpers:

- **`str_utils.py`** - String processing utilities
- **`mapper.py`** - Data mapping and transformation functions
- **`loom_mappings.py`** - Loom-specific data mappings
- **`test_ontologies.py`** - Ontology testing utilities

### `/run_lh.sh`
SLURM job script for running evaluations on the cluster:
- Configures SLURM parameters (memory, time, GPU allocation)
- Sets up environment variables and dependencies
- Executes NER evaluation with LLM validation

## Key Features

### Entity Recognition
- **Multi-type NER**: Extracts 8 different biomedical entity types
- **High-performance models**: Uses fine-tuned PubMedBERT models
- **Entity grounding**: Maps entities to standardized vocabularies
- **LLM validation**: Uses GPT-4o-mini for validation of uncertain extractions

### Knowledge Graph Construction
- **Multi-source integration**: Combines data from 10+ biomedical databases
- **BioCypher framework**: Uses BioCypher for standardized graph construction
- **Scalable processing**: Handles large-scale data with batch processing
- **Relationship extraction**: Identifies and summarizes entity relationships

### Evaluation Framework
- **Comprehensive evaluation**: Tests on multiple benchmark datasets
- **LLM-based validation**: Uses large language models for quality assessment
- **Performance metrics**: Detailed precision, recall, and F1-score reporting
- **Database integration**: Connects to Neo4j for real-time evaluation

## Usage

### Building the Knowledge Graph
```bash
python scripts/build_kg.py
```

### Running NER Evaluation
```bash
python scripts/evaluate_ner.py
```

### Using the NER System
```python
from information_extraction.NER import BiomedicalNER

ner = BiomedicalNER()
results = ner.extract_entities(
    text="BRCA1 mutations were studied in MCF-7 cells",
    entity_types=['gene', 'disease'],
    ground_entities=True,
    evaluate_with_llm=True
)
```

### Running on SLURM Cluster
```bash
sbatch run_lh.sh
```

## Dependencies

- **Core**: Python 3.9+, PyTorch, Transformers
- **NER**: seqeval, datasets, openai
- **Graph**: biocypher, neo4j
- **Data Processing**: pandas, pubmed_parser, nltk
- **Evaluation**: scikit-learn, tqdm

## Configuration

The system requires several configuration files and data paths:
- Graph schema configuration: `/data/graph_schema/glkb_schema_config.yaml`
- BioCypher configuration: `/data/graph_schema/glkb_biocypher_config.yaml`
- Data directories for various biomedical databases
- OpenAI API key for LLM-based evaluation

## Data Sources

The system integrates data from:
- PubMed (articles, abstracts, metadata)
- dbSNP (genetic variants)
- Reactome (pathways)
- Gene Ontology (annotations)
- PrimeKG (knowledge graph)
- GWAS Catalog (genetic associations)
- Various biomedical vocabularies

## Performance

- **NER Models**: State-of-the-art BERT-based models fine-tuned on biomedical text
- **Evaluation**: Comprehensive testing on multiple benchmark datasets
- **Scalability**: Designed for processing large-scale biomedical literature
- **Accuracy**: LLM-based validation ensures high-quality extractions
