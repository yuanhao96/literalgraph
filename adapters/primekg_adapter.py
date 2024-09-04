# using https://github.com/biocypher/project-template/blob/main/template_package/adapters/example_adapter.py as blueprint

import string
from enum import Enum, auto
from itertools import chain
from biocypher._logger import logger
from adapters import Adapter, Node, Edge
from utils.str_utils import escape_text
from utils.mapper import biomart_mapper, drugbank_mapper
from entity_mapping.gilda_grounders import Disease_Grounder, Chemical_Grounder
import pandas as pd
logger.debug(f"Loading module {__name__}.")

LABEL_MAPPING = {
    'protein_protein': 'gene_gene', 
    'drug_protein': 'chemical_gene', 
    'contraindication': 'chemical_disease',
    'indication': 'chemical_disease',
    'off-label use': 'chemical_disease',
    'drug_drug': 'chemical_chemical', 
    'phenotype_protein': 'gene_disease',
    'disease_phenotype_negative': 'disease_disease',
    'disease_phenotype_positive': 'disease_disease', 
    'disease_protein': 'gene_disease', 
    'drug_effect': 'chemical_disease',
    'molfunc_protein': 'gene_go', 
    'cellcomp_protein': 'gene_go', 
    'bioprocess_protein': 'gene_go', 
    'exposure_protein': 'chemical_gene', 
    'exposure_disease': 'chemical_disease',
    'exposure_bioprocess': 'chemical_go', 
    'exposure_molfunc': 'chemical_go', 
    'exposure_cellcomp': 'chemical_go', 
    'pathway_protein': 'gene_to_pathway_association',
    'anatomy_protein_present': 'gene_anatomy',
    'anatomy_protein_absent': 'gene_anatomy'
}

GENE_MAPPER = biomart_mapper()
CHEM_MAPPER = drugbank_mapper()
DISEASE_GROUNDER = Disease_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/disease.json')
CHEM_GROUNDER = Chemical_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/chemical.json')

def ground_primekg(source, id, name):
    if source=='NCBI':
        m = GENE_MAPPER.get(str(id), 'entrez')
        if m:
            return f'hgnc:{m}'
    elif source=='DrugBank':
        m = CHEM_MAPPER.get(str(id))
        if m:
            return f'chebi:{m}'
    elif source=='HPO':
        return f'hp:{id.zfill(7)}'
    elif source in ['MONDO', 'GO', 'UBERON']:
        return f'{source.lower()}:{id.zfill(7)}'
    elif source=='REACTOME':
        return f'{source.lower()}:{id}'
    elif source=='MONDO_grouped':
        terms = DISEASE_GROUNDER.ground(name)
        ms =  [m.term.get_curie() for m in terms]
        if len(ms)>0:
            return ms[0]
    elif source=='CTD':
        terms = CHEM_GROUNDER.ground(name)
        ms = [m.term.get_curie() for m in terms]
        if len(ms)>0:
            return ms[0]

class PrimeKGAdapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    pass

class PrimeKGAdapter_EdgeType(Enum):
    """
    Enum for the types of the adapter.
    """
    pass

class PrimeKGAdapter(Adapter):
    """
    PrimeKG BioCypher adapter. Generates nodes and edges for creating a
    knowledge graph.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
    """
    def __init__(self,
        node_types: str = None,
        node_fields: str = None,
        edge_types: str = None,
        edge_fields: str = None
        ):
        self._set_types_and_fields(
            node_types, node_fields, edge_types, edge_fields
        )
        self.nodes = None
        self.edges = None
        self.data = None

    def _set_types_and_fields(
        self, node_types, node_fields, edge_types, edge_fields
    ):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type.value for type in PrimeKGAdapter_NodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field.value
                for field in chain(
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type.value for type in PrimeKGAdapter_EdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field.value 
                for field in chain(
                )
            ]

    def get_nodes(self, file:str = None):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """
        logger.info("Generating nodes.")
        if file:
            self.load_data(file=file)
        elif not self.data:
            raise Exception('Please provide a PrimeKG file, or run load_data first!')
        if not self.nodes:
            self.nodes = []

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())
    
    def get_edges(self, file:str = None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.
        """
        logger.info("Generating edges.")
        if file:
            self.load_data(file=file)
        elif not self.data:
            raise Exception('Please provide a PrimeKG file, or run load_data first!')
        if not self.edges:
            self.edges = []
        
        for d in self.data:
            self.edges.append(
                Edge(
                source=d.get('head'),
                target=d.get('tail'),
                label=d.get('label'),
                properties=d
                ))
        
        for edge in self.edges:
            yield (edge.get_id(), edge.get_source(), edge.get_target(), edge.get_label(), edge.get_properties())

    def load_data(self, data:str):
        """
        Parse processed PrimeKG
        """
        logger.info("Loading PrimeKG from disk.")
        self.data = []

        # data
        df = pd.read_csv(data).astype(str)
        df = df[df['display_relation']!='parent-child']
        df.columns = ['label', 'type', 'x_index', 'x_id', 'x_type', 'x_name', 'x_source', 'y_index', 'y_id', 'y_type', 'y_name', 'y_source']
        df['head'] = df.apply(lambda x: ground_primekg(x['x_source'], x['x_id'], x['x_name']), axis=1)
        df['tail'] = df.apply(lambda x: ground_primekg(x['y_source'], x['y_id'], x['y_name']), axis=1)
        df = df[['head', 'tail', 'label', 'type']]
        df['label'] = df['label'].apply(LABEL_MAPPING.get)
        df['source'] = 'primekg'
        self.data += list(df.drop_duplicates().dropna().to_dict(orient='index').values())
        return self