# using https://github.com/biocypher/project-template/blob/main/template_package/adapters/example_adapter.py as blueprint

import pyobo
import string
from enum import Enum, auto
from itertools import chain
from biocypher._logger import logger
from adapters import Adapter, Node, Edge
from utils.str_utils import escape_text
import pandas as pd

logger.debug(f"Loading module {__name__}.")
PREFIXES = {
    'gene': ["hgnc"],
    'disease': ['mondo', 'doid'],
    'phenotype': ["hp"],
    'chemical': ["chebi"],
    'anatomy': ['bto', 'uberon', 'efo', 'cl'],
    'organism': ["ncbitaxon"],
    'mesh_term': ['mesh']
}

def get_om_curie(id):
    target_prefixes = set([j for i in PREFIXES.values() for j in i])
    lst = id.split('_')
    if len(lst)==2:
        if lst[0]=='hugo.owl#hgnc': # special rules for hgnc
            lst[0] = 'hgnc'
        if lst[0].lower() in target_prefixes:
            return f"{lst[0].lower()}:{lst[1]}"

class OntologyAdapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    GENE = "gene"
    DISEASE = "disease"
    CHEMICAL = "chemical"
    ANATOMY = "anatomy"
    ORGANISM = "organism"
    PHENOTYPE = "phenotype"
    MESH = "mesh_term"

class OntologyAdapter_Field(Enum):
    """
    Define possible fields the adapter can provide for all nodes.
    """

    NAME = "name"
    SYNONYMS = "synonyms"
    DESCRIPTION = "description"
    SOURCE = "source"

class OntologyAdapter_EdgeType(Enum):
    """
    Enum for the types of the adapter.
    """
    HIER = "hierarchical_structure"

class OntologyAdapter_Hier_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for hierarchical_structure edges.
    """
    SOURCE = "source"
    TYPE = "type"

class OMAdapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    pass

class OMAdapter_EdgeType(Enum):
    """
    Enum for the types of the adapter.
    """
    ONTOLOGY_MAPPING = "ontology_mapping"

class OMAdapter_OM_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for ontology_mapping edges.
    """
    SOURCE = "source"
    SCORE = "score"

class OntologyAdapter(Adapter):
    """
    OBO Ontology BioCypher adapter. Generates nodes and edges for creating a
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
            self.node_types = [type.value for type in OntologyAdapter_NodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field.value
                for field in chain(
                    OntologyAdapter_Field,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type.value for type in OntologyAdapter_EdgeType]

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
            raise Exception('Please provide a BERN2 annotation, or run load_data first!')
        if not self.nodes:
            self.nodes = []
        
        for d in self.data['nodes']:
            self.nodes.append(OBOConcept(
                id = d["curie"],
                label=d['label'],
                fields=self.node_fields,
                properties=d
            ))
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
            raise Exception('Please provide a BERN2 annotation, or run load_data first!')
        if not self.edges:
            self.edges = []

        for d in self.data['edges']:
            self.edges.append(Hier(
                source=d.get('head'),
                target=d.get('tail'),
                properties=d
            ))
        for edge in self.edges:
            yield (edge.get_id(), edge.get_source(), edge.get_target(), edge.get_label(), edge.get_properties())

    def load_data(self, prefiexes: dict=None):
        """
        Parse obo ontology
        """
        logger.info("Parsing obo onotlogies")
        if not prefiexes:
             prefiexes = PREFIXES
        self.data = {
            'nodes': [],
            'edges': []
        }
        

        for n in self.node_types:
            if prefiexes.get(n):
                for prefix in prefiexes.get(n):
                    onto = pyobo.get_ontology(prefix)
                    # mapping = pyobo.get_id_name_mapping(prefix)
                    mapping = onto.get_id_name_mapping()
                    # synonyms = pyobo.get_id_synonyms_mapping(prefix)\
                    synonyms = onto.get_id_synonyms_mapping()
                    # descriptions = pyobo.get_id_definition_mapping(prefix)
                    descriptions = onto.get_id_definition_mapping()
                    df = onto.get_relations_df()
                    
                    for identifier, name in mapping.items():
                        curie = f'{prefix}:{identifier}'
                        dat = {
                            'curie': curie,
                            'label': n,
                            'name': name,
                            'synonyms': synonyms.get(identifier),
                            'description': descriptions.get(identifier)
                                          }
                        self.data['nodes'].append(dat)
                    
                    df.columns = ['head', 'rel_ns', 'type', 'source', 'tail']
                    df = df[df['source']==prefix]
                    df['head'] = df['head'].apply(lambda x: f'{prefix}:{x}')
                    df['tail'] = df['tail'].apply(lambda x: f'{prefix}:{x}')
                    self.data['edges'] += list(df.to_dict(orient='index').values())

        return self

class OMAdapter(Adapter):
    """
    Ontology mapping BioCypher adapter. Generates nodes and edges for creating a
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
            self.node_types = [type.value for type in OMAdapter_NodeType]

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
            self.edge_types = [type.value for type in OMAdapter_EdgeType]

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
                OM(
                source=d.get('head'),
                target=d.get('tail'),
                properties=d
                ))
        
        for edge in self.edges:
            yield (edge.get_id(), edge.get_source(), edge.get_target(), edge.get_label(), edge.get_properties())

    def load_data(self, data:str):
        """
        Parse processed Ontology Mapping
        """
        logger.info("Loading OM from disk.")
        self.data = []

        # data
        df = pd.read_csv(data)
        df = pd.read_csv('/nfs/turbo/umms-drjieliu/usr/xinyubao/umls_matching/database/mappings_without_dup.csv')
        df.columns = ['head', 'tail', 'type', 'score', 'source']
        df['head'] = df['head'].apply(get_om_curie)
        df['tail'] = df['tail'].apply(get_om_curie)
        df = df[['head', 'tail', 'score', 'source']]
        self.data += list(df.drop_duplicates().dropna().to_dict(orient='index').values())
        return self

class OBOConcept(Node):
    """
    obo concept nodes
    """
    def __init__(self, label:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.label = label
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is not None:
            return [i.value for i in OntologyAdapter_Field]
        else:
            return fields

    def _generate_properties(self, properties):
        prop_dict = {}
        for field in self.fields:
            f = properties.get(field)
            if f:
                if isinstance(f, str):
                    prop_dict[field] = escape_text(f)
                elif isinstance(f, list):
                    prop_dict[field] = [escape_text(t) for t in f]
                else:
                    prop_dict[field] = properties.get(field)
        
        return prop_dict

class Hier(Edge):
    """
    hierarchical_structure edges
    """
    def __init__(self, source:str, target:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = "hierarchical_structure"
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in OntologyAdapter_Hier_EdgeField]
        else:
            return fields

    def _generate_properties(self, properties):
        prop_dict = {}
        for field in self.fields:
            f = properties.get(field)
            if f:
                if isinstance(f, str):
                    prop_dict[field] = escape_text(f)
                else:
                    prop_dict[field] = properties.get(field)
        
        return prop_dict

class OM(Edge):
    """
    ontology_mapping edges
    """
    def __init__(self, source:str, target:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = "ontology_mapping"
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in OMAdapter_OM_EdgeField]
        else:
            return fields

    def _generate_properties(self, properties):
        prop_dict = {}
        for field in self.fields:
            f = properties.get(field)
            if f:
                if isinstance(f, str):
                    prop_dict[field] = escape_text(f)
                else:
                    prop_dict[field] = properties.get(field)
        
        return prop_dict