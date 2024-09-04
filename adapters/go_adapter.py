# using https://github.com/biocypher/project-template/blob/main/template_package/adapters/example_adapter.py as blueprint

import pyobo
import string
from enum import Enum, auto
from itertools import chain
from biocypher._logger import logger
from adapters import Adapter, Node, Edge
from utils.str_utils import escape_text
import pandas as pd
from collections import defaultdict

logger.debug(f"Loading module {__name__}.")
ID = {
    'biological_process': '0008150',
    'cellular_component': '0005575',
    "molecular_function": '0003674'
}

class GOAdapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    CC = "cellular_component" # GO:0005575
    BP = "biological_process" # GO:0008150
    MF = "molecular_function" # GO:0003674

class GOAdapter_Field(Enum):
    """
    Define possible fields the adapter can provide for all nodes.
    """

    NAME = "name"
    DESCRIPTION = "description"
    SOURCE = "source"

class GOAdapter_EdgeType(Enum):
    """
    Enum for the types of the adapter.
    """
    HIER = "hierarchical_structure"

class GOAdapter_Hier_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for hierarchical_structure edges.
    """
    SOURCE = "source"
    TYPE = "type"

class GOAdapter(Adapter):
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
            self.node_types = [type.value for type in GOAdapter_NodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field.value
                for field in chain(
                    GOAdapter_Field,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type.value for type in GOAdapter_EdgeType]

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
        
        for c in self.node_types:
            for d in self.data[c]:
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

    def load_data(self, ids: dict=None):
        """
        Parse obo ontology
        """
        logger.info("Parsing obo onotlogies")
        if not ids:
            ids = ID
        self.data = defaultdict(list)

        onto = pyobo.get_ontology('go')
        mapping = onto.get_id_name_mapping()
        synonyms = onto.get_id_synonyms_mapping()
        descriptions = onto.get_id_definition_mapping()
        df = onto.get_relations_df()
        
        for n in self.node_types:
            for identifier in onto.descendants(ids[n]):
                curie = f'go:{identifier}'
                dat = {
                    'curie': curie,
                    'label': n,
                    'name': mapping.get(identifier),
                    'synonyms': synonyms.get(identifier),
                    'description': descriptions.get(identifier)
                                    }
                self.data[n].append(dat)
        
        df.columns = ['head', 'rel_ns', 'type', 'source', 'tail']
        df = df[df['source']=='go']
        df['head'] = df['head'].apply(lambda x: f'go:{x}')
        df['tail'] = df['tail'].apply(lambda x: f'go:{x}')
        self.data['edges'] += list(df.to_dict(orient='index').values())
                    

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
            return [i.value for i in GOAdapter_Field]
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
            return [i.value for i in GOAdapter_Hier_EdgeField]
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