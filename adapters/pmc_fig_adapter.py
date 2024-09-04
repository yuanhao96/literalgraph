# data: /nfs/turbo/umms-drjieliu/proj/medlineKG/data/figure_json_by_article/pmcimage_data.json

# using https://github.com/biocypher/project-template/blob/main/template_package/adapters/example_adapter.py as blueprint

import string
from enum import Enum, auto
from itertools import chain
from biocypher._logger import logger
import re
import requests
from adapters import Adapter, Node, Edge
import json
from utils.str_utils import escape_text

logger.debug(f"Loading module {__name__}.")

class PMCFigAdapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    FIGURE = "publication_figure"


class PMCFigAdapter_Figure_Field(Enum):
    """
    Define possible fields the adapter can provide for figures.
    """
    URL = "url"
    CAPTION = "caption"
    MENTION = "mention"

class PMCFigAdapter_EdgeType(Enum):
    """
    Enum for the types of the adapter.
    """
    CONTAIN_FIGURE = "contain_figure"

class PMCFigAdapter_ContainFigure_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for contain_figure edges.
    """
    SOURCE = "source"

class PMCFigAdapter(Adapter):
    """
    PMC Figure BioCypher adapter. Generates nodes and edges for creating a
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
            self.node_types = [type.value for type in PMCFigAdapter_NodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field.value
                for field in chain(
                    PMCFigAdapter_Figure_Field,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type.value for type in PMCFigAdapter_EdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field.value 
                for field in chain(
                    PMCFigAdapter_ContainFigure_EdgeField,
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
        
        for data in self.data:
            self.nodes.append(PMCFigure(
                id = data['fid'],
                fields=self.node_fields,
                properties=data
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
            raise Exception('Please provide a OPENI dataset, or run load_data first!')
        if not self.edges:
            self.edges = []

        for data in self.data:
            self.edges.append(ContainFigure(
                source=data['article'],
                target=data['fid'],
                properties={"source": "OPENI"}
                ))
            
        for edge in self.edges:
            yield (edge.get_id(), edge.get_source(), edge.get_target(), edge.get_label(), edge.get_properties())

    def load_data(self, file:str):
        """
        Parse OPENI PMC figure annotation
        """
        logger.info("Loading openi dataset from disk.")
        self.data = [dat for dat in self.parse_pmc_json(open(file))]
        return self
    
    def parse_pmc_json(self, f):
        base_url = "https://openi.nlm.nih.gov"
        for dat in json.load(f):
            try:
                data = dict()
                data['article'] = f"pmid{dat['pmid']}"
                data['url'] = f"{base_url}{dat['imgLarge']}"
                data['fid'] = str(dat['image']['id'])
                data['source'] = 'OPENI'
            except:
                continue
            try:
                if dat['image'].get('caption'):
                    data['caption'] = str(dat['image'].get('caption'))
                res = requests.get(f"{base_url}{dat['detailedQueryURL']}")
                if res.ok:
                    details = res.json()
                    details = details['list'][0]
                    if details['image'].get('mention'):
                        data['mention'] = str(details['image'].get('mention'))
            except:
                yield data
                        
            yield data

class PMCFigure(Node):
    """
    PMCFigure nodes
    """
    def __init__(self, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.label = 'publication_figure'
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in PMCFigAdapter_Figure_Field]
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

class ContainFigure(Edge):
    """
    contain_figure edges
    """
    def __init__(self, source:str, target:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = "contain_figure"
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in PMCFigAdapter_ContainFigure_EdgeField]
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