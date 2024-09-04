# using https://github.com/biocypher/project-template/blob/main/template_package/adapters/example_adapter.py as blueprint

import string
from enum import Enum, auto
from itertools import chain
from biocypher._logger import logger
from adapters import Adapter, Node, Edge
import re
from entity_mapping.gilda_grounders import Gene_Grounder
from utils.str_utils import escape_text
from collections import defaultdict

logger.debug(f"Loading module {__name__}.")

LOCAL_GROUNDER = True
if LOCAL_GROUNDER:
    GROUNDER = Gene_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/gene.json')
else:
    GROUNDER = Gene_Grounder()

class EventAdapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    MENTION = "genomic_mention"


class EventAdapter_Event_Field(Enum):
    """
    Define possible fields the adapter can provide for events.
    """
    TEXT = "text"
    CHAR_START = "char_start"
    CHAR_END = "char_end"
    SOURCE = "source"

class EventAdapter_Mention_Field(Enum):
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

class EventAdapter_EdgeType(Enum):
    """
    Enum for the types of the adapter.
    """
    IN_EVENT = "in_event"

class EventAdapter_InEvent_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for contain_mention edges.
    """
    TEXT = "text"
    SOURCE = "source"
    CHAR_START = "char_start"
    CHAR_END = "char_end"
    ROLE = "role"
    FACTUALITY = "factuality"

class EventAdapter_ContainMention_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for contain_mention edges.
    """
    SOURCE = "source"

class EventAdapter_MapTo_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for mapped_to edges.
    """
    SOURCE = "source"

class EventAdapter_ContainTerm_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for contain_term edges.
    """
    TEXT = "text"
    TYPE = "type"
    SOURCE = "source"

class EventAdapter(Adapter):
    """
    Genomic event BioCypher adapter. Generates nodes and edges for creating a
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
            self.node_types = [type.value for type in EventAdapter_NodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field.value
                for field in chain(
                    EventAdapter_Event_Field,
                    EventAdapter_Mention_Field,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type.value for type in EventAdapter_EdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field.value 
                for field in chain(
                    EventAdapter_InEvent_EdgeField,
                    EventAdapter_ContainTerm_EdgeField,
                    EventAdapter_ContainMention_EdgeField,
                    EventAdapter_MapTo_EdgeField,
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
        elif self.data is None:
            raise Exception('Please provide a DEM annotation, or run load_data first!')
        if not self.nodes:
            self.nodes = []
        for dat in self.data:
            for data in dat['event'].values():
                self.nodes.append(Event(
                    id = data['id'],
                    label=data['label'],
                    fields=self.node_fields,
                    properties=data
                ))
            for data in dat['term'].values():
                self.nodes.append(GenomicMention(
                    id = data['id'],
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
        elif self.data is None:
            raise Exception('Please provide a BERN2 annotation, or run load_data first!')
        if not self.edges:
            self.edges = []

        for dat in self.data:
            for i in dat['links'].values():
                for data in i:
                    self.edges.append(InEvent(
                        source=data['head'],
                        target=data['tail'],
                        properties=data
                    ))
            for i in dat['mentions'].values():
                for data in i:
                    self.edges.append(ContainMention(
                        source=data['head'],
                        target=data['tail'],
                        properties=data
                    ))
            for i in dat['mapto'].values():
                for data in i:
                    self.edges.append(MapTo(
                        source=data['head'],
                        target=data['tail'],
                        properties=data
                    ))
            for i in dat['containterm'].values():
                for data in i:
                    self.edges.append(ContainTerm(
                        source=data['head'],
                        target=data['tail'],
                        properties=data
                    ))

        for edge in self.edges:
            yield (edge.get_id(), edge.get_source(), edge.get_target(), edge.get_label(), edge.get_properties())

    def load_data(self, file:str):
        """
        Parse DEM annotation
        """
        logger.info("Loading DEM annotation from disk.")
        self.data = self.parse_dem_file(open(file)).values()
        return self
    
    def ground_mention(self, text):
        """
        map genomic mentions to glkb vocabulary
        """
        if GROUNDER:
            terms = GROUNDER.ground(text)
            return [m.term.get_curie() for m in terms]
        return []
    
    def process_text(self, text):
        text = re.sub(r' ([^a-zA-Z\d\s]+?)', '\g<1>', text)
        text = re.sub(r'([^a-zA-Z\d\s]+?) ', '\g<1>', text)
        return text

    def parse_dem_file(self, f):
        event_types = [
            "Gene_expression",
            "Transcription",
            "Protein_catabolism",
            "Phosphorylation",
            "Localization",
            "Binding",
            "Regulation",
            "Positive_regulation",
            "Negative_regulation"
        ]

        # initialize
        event2event_id = {}
        meta_event = set()
        all_data = {}
        pmid = None

        for l in f:
            if l.startswith('###'):
                if pmid: # not first line
                    if len(data['links'])>0: # have event
                        all_data[pmid] = data
                pmid = l.strip()[3:]
                data = {
                    "term": {},
                    "event": {},
                    "links": defaultdict(list),
                    "mentions": defaultdict(list),
                    "mapto": defaultdict(list),
                    "containterm": defaultdict(list)
                }
                meta_event = set()
            elif len(l.strip())>0:
                lst = l.strip().split('\t')
                if l.startswith('T'): # entity
                    lst[4] = self.process_text(lst[4])
                    if lst[1] in event_types:
                        data['event'][lst[0]] = {
                            "id": f"dem_event_{pmid}_{lst[0]}",
                            "label": f"event_{lst[1].lower()}",
                            "text": lst[4],
                            "source": "DEM",
                        }
                        data["mentions"][lst[0]].append({
                            "head": f"pmid{pmid}",
                            "tail": f"dem_event_{pmid}_{lst[0]}",
                            "source": "DEM"
                        })
                    else:
                        data['term'][lst[0]] = {
                            "id": f"dem_ent_{pmid}_{lst[0]}",
                            "text": lst[4],
                            "map": self.ground_mention(lst[4]),
                            "source": "DEM"
                        }
                        data["mentions"][lst[0]].append({
                            "head": f"pmid{pmid}",
                            "tail": f"dem_ent_{pmid}_{lst[0]}",
                            "source": "DEM"
                        })
                        for m in data['term'][lst[0]]['map']:
                            data["mapto"][lst[0]].append({
                                "head": f"dem_ent_{pmid}_{lst[0]}",
                                "tail": m,
                                "source": "DEM"
                            })
                            data["containterm"][lst[0]].append({
                                "head": f"pmid{pmid}",
                                "tail": m,
                                "source": "DEM"
                            })
                elif l.startswith('E'): # event
                    lst = lst[0:1] + lst[1].split(' ')
                    event_id = lst[1].split(':')[1]
                    event = data['event'].get(event_id)
                    event2event_id[lst[0]] = event_id
                    if event:
                        for r in lst[2:]:
                            role, term_id = r.split(':')
                            if term_id.startswith('T'): # primary event
                                term = data['term'].get(term_id)
                                if term:
                                    if event.get('terms'):
                                        event['terms'].add(term_id)
                                    else:
                                        event['terms'] = set([term_id])
                                    pair = (event_id, term_id)
                                    if pair not in meta_event: # new event rel
                                        # for m in term['map']:
                                        #     data["links"][lst[0]].append({
                                        #         "head": m,
                                        #         "tail": event['id'],
                                        #         "source": "DEM",
                                        #         "role": role,
                                        #         "text": term['text']
                                        #     })
                                        data["links"][lst[0]].append({
                                            "head": term['id'],
                                            "tail": event['id'],
                                            "source": "DEM",
                                            "role": role,
                                            "text": term['text']
                                        })
                                    meta_event.add(pair)
                            elif term_id.startswith('E'): # secondary event
                                # term = data['event'].get(event2event_id[term_id])
                                # if term.get('terms'):
                                #     for _id in term.get('terms'):
                                #         term = data['term'].get(_id)
                                #         if term:
                                #             if event.get('terms'):
                                #                 event['terms'].add(_id)
                                #             else:
                                #                 event['terms'] = set([_id])
                                #             pair = (event_id, term_id)
                                #             if pair not in meta_event: # new event rel
                                #                 for m in term['map']:
                                #                     data["links"][lst[0]].append({
                                #                         "head": m,
                                #                         "tail": event['id'],
                                #                         "source": "DEM",
                                #                         "role": role,
                                #                         "text": term['text']
                                #                     })
                                #             meta_event.add(pair)
                                term = data['event'].get(term_id)
                                if term:
                                    pair = (event_id, term_id)
                                    if pair not in meta_event: # new event rel
                                        data["links"][lst[0]].append({
                                            "head": term['id'],
                                            "tail": event['id'],
                                            "source": "DEM",
                                            "role": role,
                                            "text": term['text']
                                        })
                                    meta_event.add(pair)
                elif l.startswith('M'): # factuality
                    lst = lst[0:1] + lst[1].split(' ')
                    for link in data["links"][lst[2]]:
                        link['factuality'] = lst[1]
            else:
                continue
        return all_data

class Event(Node):
    """
    genomic event nodes
    """
    def __init__(self, label, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.label = label
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in EventAdapter_Event_Field]
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
            return [i.value for i in EventAdapter_Mention_Field]
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

class InEvent(Edge):
    """
    in_event edges
    """
    def __init__(self, source:str, target:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = "in_event"
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in EventAdapter_InEvent_EdgeField]
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
            return [i.value for i in EventAdapter_ContainMention_EdgeField]
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
            return [i.value for i in EventAdapter_MapTo_EdgeField]
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
            return [i.value for i in EventAdapter_ContainTerm_EdgeField]
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
            return [i.value for i in EventAdapter_ContainMention_EdgeField]
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
            return [i.value for i in EventAdapter_MapTo_EdgeField]
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