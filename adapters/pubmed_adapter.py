# using https://github.com/biocypher/project-template/blob/main/template_package/adapters/example_adapter.py as blueprint

import string
from enum import Enum, auto
from itertools import chain
from biocypher._logger import logger
from adapters import Adapter, Node, Edge
import gzip
import pubmed_parser as pp
from utils.str_utils import escape_text
from nltk.tokenize import sent_tokenize

logger.debug(f"Loading module {__name__}.")

class PubmedAdapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    ARTICLE = "pubmed_article"
    SENT = "sentence"


class PubmedAdapter_Article_Field(Enum):
    """
    Define possible fields the adapter can provide for articles.
    """

    PMID = "pmid"
    TITLE = "title"
    ABSTRACT = "abstract"
    PUBDATE = "pubdate"
    AUTHORS = "authors"
    AUTHOR_AFFILIATIONS = "author_affiliations"
    JOURNAL = "journal"
    SOURCE = "source"
    PUBTYPE = "publication_type"
    PMCID = "pmcid"
    DOI = "doi"
    PUBMEDID = "pubmedid"

class PubmedAdapter_Sentence_Field(Enum):
    """
    Define possible fields the adapter can provide for sentences.
    """
    SENTID = 'sentid'
    TEXT = 'text'

class PubmedAdapter_EdgeType(Enum):
    """
    Enum for the types of the pubmed adapter.
    """

    PUBLISHED_IN = "published_in"
    CITE = "cite"
    CONTAIN_TERM = "contain_term"
    CONTAIN_SENT = "contain_sentence"


class PubmedAdapter_PublishedIn_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for published_in edges.
    """
    pass

class PubmedAdapter_Cite_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for cite edges.
    """
    pass

class PubmedAdapter_ContainTerm_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for topic edges.
    """
    SOURCE = "source"

class PubmedAdapter_ContainSent_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for sentence edges.
    """
    pass

class PubmedAdapter(Adapter):
    """
    PubMed BioCypher adapter. Generates nodes and edges for creating a
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
        self.dicts = None
        
    
    def _set_types_and_fields(
        self, node_types, node_fields, edge_types, edge_fields
    ):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type.value for type in PubmedAdapter_NodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field.value
                for field in chain(
                    PubmedAdapter_Article_Field,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type.value for type in PubmedAdapter_EdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field.value 
                for field in chain(
                    PubmedAdapter_PublishedIn_EdgeField,
                    PubmedAdapter_Cite_EdgeField,
                    PubmedAdapter_ContainTerm_EdgeField
                )
            ]

    def get_nodes(self, pubmed_xml:str = None):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """
        logger.info("Generating nodes.")
        if pubmed_xml:
            self.load_data(file=pubmed_xml)
        elif not self.dicts:
            raise Exception('Please provide a pubmed xml, or run load_data first!')
        if not self.nodes:
            self.nodes = []

        for article in self.dicts:
            article_info = self.article_node(article)
            article_info['pubtype'] = self.get_pubtype(article)
            self.nodes.append(PubmedArticle(
                id = f"pmid{article_info['pmid']}",
                fields=self.node_fields,
                properties=article_info
            ))

            for sent_info in self.sentence_node(article):
                self.nodes.append(Sentence(
                id = sent_info.get('sentid'),
                fields=self.node_fields,
                properties=sent_info
                ))
        
        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())
    
    def get_edges(self, pubmed_xml:str = None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.
        """
        logger.info("Generating edges.")
        if pubmed_xml:
            self.load_data(file=pubmed_xml)
        elif not self.dicts:
            raise Exception('Please provide a pubmed xml, or run load_data first!')
        if not self.edges:
            self.edges = []
        
        # fix this, connect with edge field classes above
        for article in self.dicts:
            pmid = f"pmid{article['pmid']}"
            for i, topic in enumerate(self.get_mesh(article)):
                self.edges.append(Edge(
                    # id=f"Pmid{pmid}2Topic{i}",
                    source=pmid,
                    target=f"mesh:{topic}",
                    label=PubmedAdapter_EdgeType.CONTAIN_TERM.value,
                    properties={'source': "PubMed"}
                ))
            if self.get_references(article):
                for i, citation in enumerate(self.get_references(article)):
                    if citation.isnumeric(): # valid pmid
                        self.edges.append(Edge(
                            # id=f"Pmid{pmid}2Ref{i}",
                            source=pmid,
                            target=f"pmid{citation}",
                            label=PubmedAdapter_EdgeType.CITE.value,
                            properties={'source': "PubMed"}
                        ))
            if self.get_journal_id(article):
                self.edges.append(Edge(
                        # id=f"Pmid{pmid}2Journal",
                        source=pmid,
                        target=f"nlmid{self.get_journal_id(article)}",
                        label=PubmedAdapter_EdgeType.PUBLISHED_IN.value,
                        properties={}
                    ))
            for sent_info in self.sentence_node(article):
                self.edges.append(Edge(
                        source=pmid,
                        target=sent_info.get('sentid'),
                        label=PubmedAdapter_EdgeType.CONTAIN_SENT.value,
                        properties={}
                    ))
            
        for edge in self.edges:
            yield (edge.get_id(), edge.get_source(), edge.get_target(), edge.get_label(), edge.get_properties())

    def get_pmids(self, pubmed_xml:str = None):
        """
        Returns a list of str(PMID) inside the xml file.
        """
        logger.info("Listing PMIDs.")
        if pubmed_xml:
            self.load_data(file=pubmed_xml)
        elif not self.dicts:
            raise Exception('Please provide a pubmed xml, or run load_data first!')
        
        # print(self.dicts)
        pmids = [str(article['pmid']) for article in self.dicts]
        return pmids
    
    def load_data(self, file:str):
        """
        Parse PubMed primary source
        """
        logger.info("Loading PubMed data from disk.")
        # if file.endswith('gz'):
        #     file = gzip.open(file)
        # else:
        #     file = open(file)

        dicts_out = pp.parse_medline_xml(
        file,
        year_info_only=True,
        author_list=True,
        reference_list=True,) # return list of dictionary
        self.dicts = dicts_out
        return self

    def sentence_node(self, article):
        try:
            sents = []
            title = article['title']
            sents.append(title)
            _id = article['pmid']
            abstract = article['abstract']
            if isinstance(abstract, str):
                sents += sent_tokenize(abstract)
            for i, sent in enumerate(sents):
                yield {'sentid': f'pmid{_id}_{i}', 'text':sent}
        except KeyError:
            return sents

    def article_node(self, article):
        try:
            _id = article['pmid']
            doi = article['doi']
            pmc = article['pmc']
            title = article['title']
            date = article['pubdate']
            abstract = article['abstract']
            try:
                journal = article['journal'].split(' = ')[0]
            except:
                journal = 'Unknown Journal'
            try:
                authors = '|'.join([a['forename']+' '+a['lastname'] for a in article['authors']])
                author_ids = '|'.join(list(set(a['identifier'] for a in article['authors']).difference({''})))
                author_affiliations = '|'.join(list(set(a['affiliation'] for a in article['authors']).difference({''})))
            except KeyError:
                authors = ''
                author_ids = ''
            # remove escape chars in title
            escapes = ''.join([chr(char) for char in range(1, 32)])
            translator = str.maketrans('', '', escapes)
            title = title.translate(translator)
            return {'pmid':_id, 'pubmedid':_id, 'doi':doi, 'pmc':pmc, 'title':title, 'pubdate':date, 'authors':authors, 'author_ids':author_ids, 'author_affiliations':author_affiliations, 'abstract':abstract, 'journal':journal, 'source':'PubMed'}
        
        except KeyError:
            return None
    
    def get_journal(self, article):
        try:
            return article['journal']
        except KeyError:
            return None
    
    def get_journal_id(self, article):
        if article.get('nlm_unique_id'):
            return article.get('nlm_unique_id')

    def get_pubtype(self, article):
        try:
            return [p.split(':')[-1] for p in article['publication_types'].split(';')]
        except KeyError:
            return None
    
    def get_mesh(self, article):
        try:
            if len(article['mesh_terms']) > 0:
                return [m.split(':')[0] for m in article['mesh_terms'].split('; ')]
            else:
                return []
        except KeyError:
            return []
    
    def get_references(self, article):
        try:
            return [ref['pmid'] for ref in article['references']]
        except KeyError:
            return []

class PubmedArticle(Node):
    """
    PubMed article nodes
    """
    def __init__(self, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.label = 'pubmed_article'
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is not None:
            return [i.value for i in PubmedAdapter_Article_Field]
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

class Sentence(Node):
    """
    PubMed sentence nodes
    """
    def __init__(self, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.label = 'sentence'
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is not None:
            return [i.value for i in PubmedAdapter_Sentence_Field]
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


class PublishedIn(Edge):
    """
    published_in edges
    """
    def __init__(self, source:str, target:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = "published_in"
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in PubmedAdapter_PublishedIn_EdgeField]
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
            return [i.value for i in PubmedAdapter_ContainTerm_EdgeField]
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
