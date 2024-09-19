# using https://github.com/biocypher/project-template/blob/main/template_package/adapters/example_adapter.py as blueprint

import string
from enum import Enum, auto
from itertools import chain
from biocypher._logger import logger
from adapters import Adapter, Node, Edge
import re
import json
from entity_mapping.gilda_grounders import Gene_Grounder, Disease_Grounder, Chemical_Grounder, Organism_Grounder, Anatomy_Grounder
from utils.str_utils import escape_text
from nltk.tokenize import sent_tokenize

logger.debug(f"Loading module {__name__}.")

LOCAL_GROUNDER = True
if LOCAL_GROUNDER:
    GROUNDERS = {
            'gene': Gene_Grounder(prefixes=None, file='configs/custom_grounders/gene.json'),
            'DNA': Gene_Grounder(prefixes=None, file='configs/custom_grounders/gene.json'),
            'disease': Disease_Grounder(prefixes=None, file='configs/custom_grounders/disease.json'),
            'drug': Chemical_Grounder(prefixes=None, file='configs/custom_grounders/chemical.json'),
            'species': Organism_Grounder(prefixes=None, file='configs/custom_grounders/organism.json'),
            'cell_type': Anatomy_Grounder(prefixes=None, file='configs/custom_grounders/anatomy.json'),
            'cell_line': Anatomy_Grounder(prefixes=None, file='configs/custom_grounders/anatomy.json'),
        }
else:
    GROUNDERS = {
            'gene': Gene_Grounder(),
            'DNA': Gene_Grounder(),
            'disease': Disease_Grounder(),
            'drug': Chemical_Grounder(),
            'species': Organism_Grounder(),
            'cell_type': Anatomy_Grounder(),
            'cell_line': Anatomy_Grounder(),
        }

class BERN2Adapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    MENTION = "genomic_mention"


class BERN2Adapter_Mention_Field(Enum):
    """
    Define possible fields the adapter can provide for genomic mentions.
    """

    TEXT = "text"
    CHAR_START = "char_start"
    CHAR_END = "char_end"
    TYPE = "type"
    SOURCE = "source"
    NORM_NAME = "normalized_name"
    PROB = "prob"

class BERN2Adapter_EdgeType(Enum):
    """
    Enum for the types of the adapter.
    """
    CONTAIN_MENTION = "contain_mention"
    MAP_TO = "mapped_to"
    CONTAIN_TERM = "contain_term"

class BERN2Adapter_ContainMention_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for contain_mention edges.
    """
    SOURCE = "source"

class BERN2Adapter_MapTo_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for mapped_to edges.
    """
    SOURCE = "source"
    SCORE = "score"

class BERN2Adapter_ContainTerm_EdgeField(Enum):
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

class BERN2Adapter(Adapter):
    """
    BERN2 BioCypher adapter. Generates nodes and edges for creating a
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
            self.node_types = [type.value for type in BERN2Adapter_NodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field.value
                for field in chain(
                    BERN2Adapter_Mention_Field,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type.value for type in BERN2Adapter_EdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field.value 
                for field in chain(
                    BERN2Adapter_ContainTerm_EdgeField,
                    BERN2Adapter_ContainMention_EdgeField,
                    BERN2Adapter_MapTo_EdgeField,
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
            self.nodes.append(GenomicMention(
                id = data['eid'],
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
            raise Exception('Please provide a BERN2 annotation, or run load_data first!')
        if not self.edges:
            self.edges = []

        for data in self.data:
            pmid = f"pmid{data['pmid']}"
            self.edges.append(ContainMention(
                source=pmid,
                target=data['eid'],
                properties={"source": "BERN2"}
                ))
            if "sentence" in data:
                self.edges.append(ContainMention(
                    source=data['sentence'],
                    target=data['eid'],
                    properties={"source": "BERN2"}
                    ))
            
            mappings = self.ground_mention(data['text'], data['type'])
            for m, score in mappings[:1]:
                self.edges.append(ContainTerm(
                    source=pmid,
                    target=m,
                    properties=data
                ))
                self.edges.append(MapTo(
                    source=data['eid'],
                    target=m,
                    properties={"source": "BERN2",
                    "score": score
                    }
                ))

            # mutations
            if data['type'] == 'mutation':
                rsid = None
                if 'rs' in data['normalized_name'].lower():
                    pat1 = re.compile(r'RS#:([0-9]+)$')
                    pat2 = re.compile(r'rs([0-9]+)$')
                    match = re.search(pat1, data['normalized_name'])
                    if match:
                        rsid = f"rs{match.group(1)}"
                    else:
                        match = re.search(pat2, data['normalized_name'])
                        if match:
                            rsid = f"rs{match.group(1)}"
                
                if rsid:
                    self.edges.append(ContainTerm(
                        source=pmid,
                        target=rsid,
                        properties=data
                    ))
                    self.edges.append(MapTo(
                        source=data['eid'],
                        target=rsid,
                        properties={"source": "BERN2"}
                    ))

        for edge in self.edges:
            yield (edge.get_id(), edge.get_source(), edge.get_target(), edge.get_label(), edge.get_properties())

    def load_data(self, file:str):
        """
        Parse BERN2 annotation
        """
        logger.info("Loading BERN2 annotation from disk.")
        self.data = [dat for dat in self.parse_bern2_json(file)]
        return self
    
    def ground_mention(self, text, type):
        """
        map genomic mentions to glkb vocabulary
        """
        if GROUNDERS.get(type):
            terms = GROUNDERS.get(type).ground(text)
            return [(m.term.get_curie(), m.score) for m in terms]
        return []

    def parse_bern2_json(self, path):
        true = True
        false = False
        NaN = None
        null = None
        skip = []

        # Determine which sentence each span is in
        def find_sentence_for_span(span, sentence_boundaries):
            for i, (sent_start, sent_end) in enumerate(sentence_boundaries):
                if span[0] >= sent_start and span[1] <= sent_end:
                    return i
            return None  # In case the span doesn't fit in any sentence

        with open(path, "r") as f:
            data = json.load(f)
        
        for l in data:
            dic = eval(l)
            sentence_boundaries = []
            if 'annotations' not in dic.keys():
                continue  
            if 'text' in dic.keys():
                sentences = sent_tokenize(dic['text'])
                start = 0
                for sentence in sentences:
                    end = start + len(sentence)
                    sentence_boundaries.append((start, end))
                    start = end + 1  # Assuming a single character (like a space or period) between sentences
                
            for i, d in enumerate(dic['annotations']):
                data = {'source': 'BERN2'}
                data['pmid'] = str(dic['pmid'])
                data['text'] = d['mention']
                data['type'] = d['obj']
                data["char_start"] = d['span']['begin']
                data["char_end"] = d['span']['end']
                data['id'] = d['id'] # mapped id
                data['eid'] = f"ent-{data['pmid']}-{data['char_start']}-{data['char_end']}-{data['type']}"
                # find the sentences containing the spans
                sentence_index = find_sentence_for_span((d['span']['begin'], d['span']['end']), sentence_boundaries)
                if sentence_index is not None:
                    data['sentence'] = f"pmid{dic['pmid']}_{sentence_index}"
                try:
                    data['prob'] = d['prob']
                except KeyError:
                    pass
                if d['obj']=='mutation': # mutation
                    data['normalized_name'] = d['normalizedName']
                else:
                    data['normalized_name'] = d['mention'].lower()
                        
                yield data

class GenomicMention(Node):
    """
    genomic mention nodes
    """
    def __init__(self, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.label = 'genomic_mention'
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in BERN2Adapter_Mention_Field]
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
            return [i.value for i in BERN2Adapter_ContainTerm_EdgeField]
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

class ContainMention(Edge):
    """
    contain_mention edges
    """
    def __init__(self, source:str, target:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = "contain_mention"
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in BERN2Adapter_ContainMention_EdgeField]
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


class MapTo(Edge):
    """
    mapped_to edges
    """
    def __init__(self, source:str, target:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = "mapped_to"
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in BERN2Adapter_MapTo_EdgeField]
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
