# using https://github.com/biocypher/project-template/blob/main/template_package/adapters/example_adapter.py as blueprint

import string
from enum import Enum, auto
from itertools import chain
from biocypher._logger import logger
from adapters import Adapter, Node, Edge
from utils.str_utils import escape_text

logger.debug(f"Loading module {__name__}.")

class JournalAdapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    JOURNAL = "journal"


class JournalAdapter_Journal_Field(Enum):
    """
    Define possible fields the adapter can provide for journals.
    """
    TITLE = "title"
    MED_ABBR = "med_abbrevation"
    ISO_ABBR = "iso_abbrevation"
    ISSN_PRINT = "issn_print"
    ISSN_ONLINE = "issn_online"
    JRID = "jrid"


class JournalAdapter_EdgeType(Enum):
    """
    Enum for the types of the journal adapter.
    """
    pass


class JournalAdapter(Adapter):
    """
    PubMed Journal BioCypher adapter. Generates nodes and edges for creating a
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
        self.journals = None
        
    
    def _set_types_and_fields(
        self, node_types, node_fields, edge_types, edge_fields
    ):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type.value for type in JournalAdapter_NodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field.value
                for field in chain(
                    JournalAdapter_Journal_Field,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type.value for type in JournalAdapter_EdgeType]

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
        elif not self.journals:
            raise Exception('Please provide a file, or run load_data first!')
        if not self.nodes:
            self.nodes = []

        for journal in self.journals:
            self.nodes.append(Journal(
                id = f"nlmid{journal['NlmId']}",
                fields=self.node_fields,
                properties={
                    'jrid': journal.get('JrId'),
                    'med_abbrevation': journal.get('MedAbbr'),
                    'iso_abbrevation': journal.get('IsoAbbr'), 
                    'issn_print': journal.get('ISSN (Print)'),
                    'issn_online': journal.get('ISSN (Online)'),
                    'title': journal.get('JournalTitle')
                }
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
        elif not self.journals:
            raise Exception('Please provide a file, or run load_data first!')
        if not self.edges:
            self.edges = []
            
        for edge in self.edges:
            yield (edge.get_id(), edge.get_source(), edge.get_target(), edge.get_label(), edge.get_properties())

    def load_data(self, file:str):
        """
        Parse PubMed journal primary source
        """
        logger.info("Loading journal data from disk.")
        self.journals = [j for j in self.journal_reader(open(file))]
        return self
    
    def journal_reader(self, file):
        file.readline() # skip the first line
        journal = {}
        for l in file:
            if l.startswith('-'):
                journal_out = journal
                journal = {}
                yield journal_out
            else:
                lst = l.strip().split(':')
                key, val = lst[0], ''.join(lst[1:])
                journal[key.strip()] = val.strip()
    
class Journal(Node):
    """
    journal nodes
    """
    def __init__(self, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.label = 'journal'
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is not None:
            return [i.value for i in JournalAdapter_Journal_Field]
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