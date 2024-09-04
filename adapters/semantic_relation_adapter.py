# using https://github.com/biocypher/project-template/blob/main/template_package/adapters/example_adapter.py as blueprint

import string
from enum import Enum, auto
from itertools import chain
from biocypher._logger import logger
from adapters import Adapter, Node, Edge
import re
from entity_mapping.gilda_grounders import Gene_Grounder, Disease_Grounder, Chemical_Grounder, Organism_Grounder, Anatomy_Grounder
from utils.str_utils import escape_text
import pandas as pd

logger.debug(f"Loading module {__name__}.")

LOCAL_GROUNDER = True
if LOCAL_GROUNDER:
    GROUNDERS = {
            'GeneOrGeneProduct': Gene_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/gene.json'),
            'DiseaseOrPhenotypicFeature': Disease_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/disease.json'),
            'ChemicalEntity': Chemical_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/chemical.json'),
        }
else:
    GROUNDERS = {
            'GeneOrGeneProduct': Gene_Grounder(),
            'DiseaseOrPhenotypicFeature': Disease_Grounder(),
            'ChemicalEntity': Chemical_Grounder(),
        }

class SemanticRelationshipAdapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """


class SemanticRelationshipAdapter_Field(Enum):
    """
    Define possible fields the adapter can provide for genomic mentions.
    """

class SemanticRelationshipAdapter_EdgeType(Enum):
    """
    Enum for the types of the adapter.
    """
    ASSOCIATION = "associate"
    POSITIVE_CORRELATION = "positive_correlation"
    NEGATIVE_CORRELATION = "negative_correlation"
    COTREATMENT = "cotreatment"
    DRUG_INTERACTION = "drug_interaction"
    BIND = "bind"
    COMPARISON = "comparison"
    CONVERSION = "conversion"

class SemanticRelationshipAdapter_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for edges.
    """
    TEXT = "text"
    CHAR_START = "char_start"
    CHAR_END = "char_end"
    SOURCE = "source"
    TAIL_TEXT = "tail_text"
    TAIL_CHAR_START = "tail_char_start"
    TAIL_CHAR_END = "tail_char_end"
    FROM_ARTICLE = "from_article"

class SemanticRelationshipAdapter(Adapter):
    """
    BioRed relationship BioCypher adapter. Generates nodes and edges for creating a
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
            self.node_types = [type.value for type in SemanticRelationshipAdapter_NodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field.value
                for field in chain(
                    SemanticRelationshipAdapter_Field,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type.value for type in SemanticRelationshipAdapter_EdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field.value 
                for field in chain(
                    SemanticRelationshipAdapter_EdgeField,
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
            raise Exception('Please provide a semantic relationship annotation, or run load_data first!')
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
            raise Exception('Please provide a semantic relationship annotation, or run load_data first!')
        if not self.edges:
            self.edges = []

        for data in self.data:
            self.edges.append(BioRED_SemanticRelationship(
                # id=f"Pmid{pmid}2Ent{data['eid']}",
                source=data['head'],
                target=data['tail'],
                label=data['label'],
                properties=data
                ))

        for edge in self.edges:
            yield (edge.get_id(), edge.get_source(), edge.get_target(), edge.get_label(), edge.get_properties())

    def load_data(self, file:str):
        """
        Parse semantic relationship annotation
        """
        logger.info("Loading semantic relationship annotation from disk.")
        self.data = self.parse_SemanticRelationship_file(file)
        return self
    
    def ground_mention(self, text, type):
        """
        map genomic mentions to glkb vocabulary
        """
        if GROUNDERS.get(type):
            terms = GROUNDERS.get(type).ground(text)
            return [m.term.get_curie() for m in terms]
        return []

    def parse_SemanticRelationship_file(self, f):
        data = []
        df = pd.read_csv(f).dropna().drop_duplicates()
        heads = df.apply(lambda x: self.ground_mention(x['Head:string'], x['Source_type']), axis=1).to_list()
        tails = df.apply(lambda x: self.ground_mention(x['Tail:string'], x['Target_type']), axis=1).to_list()
        for hs, ts, r, ht, tt, s in zip(heads, tails, df['Type:string'], df['Head:string'], df['Tail:string'], df['Source:string']):
            for h in hs:
                for t in ts:
                    if h != t:
                        data.append({
                            'id': f"{h}-{t}-{label}-{s}",
                            'head': h,
                            'tail': t,
                            'label': r.lower(),
                            'source': 'SciBERT',
                            'from_article': f'pmid{s}',
                            'text': ht,
                            'tail_text': tt
                        })
        return data
        

class BioRED_SemanticRelationship(Edge):
    """
    contain_term edges
    """
    def __init__(self, source:str, target:str, label:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = label
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in SemanticRelationshipAdapter_EdgeField]
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
