# using https://github.com/biocypher/project-template/blob/main/template_package/adapters/example_adapter.py as blueprint

import string
from enum import Enum, auto
from itertools import chain
from biocypher._logger import logger
from adapters import Adapter, Node, Edge
from utils.str_utils import escape_text

logger.debug(f"Loading module {__name__}.")


class dbSNPAdapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    SNV = "snv"

class dbSNPAdapter_Snv_Field(Enum):
    """
    Define possible fields the adapter can provide for all nodes.
    """

    RSID = "rsid"
    REF = "ref"
    ALT = "alt"
    SOURCE = "source"

class dbSNPAdapter_EdgeType(Enum):
    """
    Enum for the types of the adapter.
    """
    pass

class dbSNPAdapter(Adapter):
    """
    dbSNP BioCypher adapter. Generates nodes and edges for creating a
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
            self.node_types = [type.value for type in dbSNPAdapter_NodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field.value
                for field in chain(
                    dbSNPAdapter_Snv_Field,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type.value for type in dbSNPAdapter_EdgeType]

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
            raise Exception('Please provide a dbSNP file, or run load_data first!')
        if not self.nodes:
            self.nodes = []
        
        for d in self.data['nodes']:
            self.nodes.append(
                SNV(
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
            raise Exception('Please provide a dbSNP file, or run load_data first!')
        if not self.edges:
            self.edges = []
        
        for edge in self.edges:
            yield (edge.get_id(), edge.get_source(), edge.get_target(), edge.get_label(), edge.get_properties())

    def load_data(self, file:str):
        """
        Parse processed dbSNP
        """
        logger.info("Loading dbSNP from disk.")
        data = {
            'nodes': [],
            'edges': []
        }
        with open(file) as f:
            for l in f:
                lst = l.strip().split('\t')
                if lst[-1] == 'True': # common snp
                    data['nodes'].append({
                        'rsid': lst[3],
                        'ref': lst[4],
                        'alt': lst[5],
                        'id': lst[3],
                        'source': 'dbSNP'
                    })
        self.data = data
        return self

class SNV(Node):
    """
    snv nodes
    """
    def __init__(self, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.label = 'snv'
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is not None:
            return [i.value for i in dbSNPAdapter_Snv_Field]
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
