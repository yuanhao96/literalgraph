# using https://github.com/biocypher/project-template/blob/main/template_package/adapters/example_adapter.py as blueprint

import string
from enum import Enum, auto
from itertools import chain
from biocypher._logger import logger
from adapters import Adapter, Node, Edge
from utils.str_utils import escape_text
from utils.mapper import biomart_mapper
import pandas as pd
logger.debug(f"Loading module {__name__}.")


class ReactomeAdapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    PATHWAY = "pathway"

class ReactomeAdapter_Pathway_Field(Enum):
    """
    Define possible fields the adapter can provide for all nodes.
    """
    NAME = 'name'
    DESC = 'description'
    SOURCE = 'source'

class ReactomeAdapter_EdgeType(Enum):
    """
    Enum for the types of the adapter.
    """
    HIER = "hierarchical_structure"
    CONTAIN_TERM = "contain_term"

class ReactomeAdapter_Hier_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for hierarchical_structure edges.
    """
    SOURCE = "source"
    TYPE = "type"

class ReactomeAdapter_G2P_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for gene_to_pathway_association edges.
    """
    SOURCE = "source"

class ReactomeAdapter_ContainTerm_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for contain_term edges.
    """
    TEXT = "text"
    CHAR_START = "char_start"
    CHAR_END = "char_end"
    TYPE = "type"
    SOURCE = "source"
    NORM_NAME = "normalized_name"
    PROB = "prob"

class ReactomeAdapter(Adapter):
    """
    Reactome BioCypher adapter. Generates nodes and edges for creating a
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
            self.node_types = [type.value for type in ReactomeAdapter_NodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field.value
                for field in chain(
                    ReactomeAdapter_Pathway_Field,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type.value for type in ReactomeAdapter_EdgeType]

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
            raise Exception('Please provide a Reactome file, or run load_data first!')
        if not self.nodes:
            self.nodes = []
        
        for d in self.data['nodes']:
            self.nodes.append(
                Pathway(
                id = d['id'],
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
            raise Exception('Please provide a Reactome file, or run load_data first!')
        if not self.edges:
            self.edges = []
        
        for d in self.data['gene2pathway']:
            self.edges.append(
                G2P(
                source=d.get('gene'),
                target=d.get('id'),
                properties=d
                ))
        
        for d in self.data['hier']:
            self.edges.append(
                Hier(
                source=d.get('head'),
                target=d.get('tail'),
                properties=d
                ))

        for d in self.data['pub2pathway']:
            self.edges.append(
                ContainTerm(
                source=d.get('pmid'),
                target=d.get('id'),
                properties=d
                ))
        
        for edge in self.edges:
            yield (edge.get_id(), edge.get_source(), edge.get_target(), edge.get_label(), edge.get_properties())

    def load_data(self, data:str, rt2gene:str, hier:str, rt2pub:str):
        """
        Parse processed Reactome
        """
        logger.info("Loading Reactome from disk.")
        self.data = {
            'nodes': [],
            'gene2pathway': [],
            'hier': [],
            'pub2pathway': []
        }
        mapper = biomart_mapper()

        # data
        df = pd.read_csv(data, sep='\t', header=None)
        df.columns = ['id', 'name', 'Species']
        df = df[df['Species']=='Homo sapiens'][['id', 'name']]
        df['id'] = df['id'].apply(lambda x: f'reactome:{x}')
        df['description'] = df['name']
        df['source'] = 'reactome'
        self.data['nodes'] += list(df.drop_duplicates().dropna().to_dict(orient='index').values())

        # pathway 2 gene
        df = pd.read_csv(rt2gene, sep='\t', header=None)
        df.columns = ['gene', 'id', 'url', 'name', 'evidence code', 'Species'] # ncbi to reactome
        df['gene'] = df['gene'].apply(lambda x: mapper.get(str(x), 'entrez'))
        df['id'] = df['id'].apply(lambda x: f'reactome:{x}')
        rt2gene = df[df['Species']=='Homo sapiens'][['gene', 'id']]
        rt2gene['source'] = 'reactome'
        self.data['gene2pathway'] += list(rt2gene.dropna().drop_duplicates().to_dict(orient='index').values())

        # hier
        conv_hier = pd.read_csv(hier, sep='\t', header=None)
        conv_hier.columns = ['tail', 'head']
        conv_hier['tail'] = conv_hier['tail'].apply(lambda x: f'reactome:{x}')
        conv_hier['head'] = conv_hier['head'].apply(lambda x: f'reactome:{x}')
        conv_hier['source'] = 'reactome'
        self.data['hier'] += list(conv_hier.dropna().drop_duplicates().to_dict(orient='index').values())

        # cite
        df = pd.read_csv(rt2pub, sep='\t', header=None)
        df.columns = ['id', 'pmid']
        df['id'] = df['id'].apply(lambda x: f'reactome:{x}')
        df['pmid'] = df['pmid'].apply(lambda x: f'pmid{x}')
        df['source'] = 'reactome'
        self.data['pub2pathway'] += list(df.dropna().drop_duplicates().to_dict(orient='index').values())

        return self

class Pathway(Node):
    """
    pathway nodes
    """
    def __init__(self, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.label = 'pathway'
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is not None:
            return [i.value for i in ReactomeAdapter_Pathway_Field]
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
            return [i.value for i in ReactomeAdapter_Hier_EdgeField]
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

class ContainTerm(Edge):
    """
    contain_term edges
    """
    def __init__(self, source:str, target:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = "contain_term"
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in ReactomeAdapter_ContainTerm_EdgeField]
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

class G2P(Edge):
    """
    gene_to_pathway_association edges
    """
    def __init__(self, source:str, target:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = "gene_to_pathway_association"
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in ReactomeAdapter_G2P_EdgeField]
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