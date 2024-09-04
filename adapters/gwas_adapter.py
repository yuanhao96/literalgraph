# using https://github.com/biocypher/project-template/blob/main/template_package/adapters/example_adapter.py as blueprint

import string
from enum import Enum, auto
from itertools import chain
from biocypher._logger import logger
from adapters import Adapter, Node, Edge
from utils.str_utils import escape_text
from entity_mapping.gilda_grounders import Gene_Grounder
from collections import defaultdict
import pandas as pd
logger.debug(f"Loading module {__name__}.")

GENE_GROUNDER = Gene_Grounder(prefixes=None, file='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/gene.json')

def ground_gene(name):
    ms = [m.term.get_curie() for m in GENE_GROUNDER.ground(name)]
    if len(ms)>0:
        return ms[0]

class GWASAdapter_NodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    pass

class GWASAdapter_EdgeType(Enum):
    """
    Enum for the types of the adapter.
    """
    pass

class GWASAdapter(Adapter):
    """
    GWAS BioCypher adapter. Generates nodes and edges for creating a
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
            self.node_types = [type.value for type in GWASAdapter_NodeType]

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
            self.edge_types = [type.value for type in GWASAdapter_EdgeType]

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
            raise Exception('Please provide a GWAS file, or run load_data first!')
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
            raise Exception('Please provide a GWAS file, or run load_data first!')
        if not self.edges:
            self.edges = []
        
        for c in self.data['edges'].values():
            for d in c:
                self.edges.append(
                    Edge(
                    source=d.get('head'),
                    target=d.get('tail'),
                    label=d.get('label'),
                    properties=d
                    ))
        
        for edge in self.edges:
            yield (edge.get_id(), edge.get_source(), edge.get_target(), edge.get_label(), edge.get_properties())

    def load_data(self, snp_gene:str, snp_trait:str):
        """
        Parse processed GWAS
        """
        logger.info("Loading GWAS from disk.")
        self.data = {
            'nodes': defaultdict(list),
            'edges': defaultdict(list)
        }

        # data
        df = pd.read_csv(snp_gene, sep='\t')
        df = df.dropna()
        df.columns = ['head', 'risk allele', 'tail']
        df['tail'] =  df['tail'].apply(ground_gene)
        df['type'] = 'SNP_intra_gene'
        df['label'] = 'variant_gene'
        df['source'] = 'gwas'
        self.data['edges']['variant_gene'] += list(df.drop_duplicates().dropna(subset=['head', 'tail']).to_dict(orient='index').values())

        def transform_trait_id(id):
            lst = id.split('/')[-1].split('_')
            return f"{lst[0].lower()}:{lst[1]}"
        df = pd.read_csv(snp_trait, sep='\t')
        df.columns = ['trait', 'tail', 'head', 'chr', 'start', 'end', 'risk allele', 'type', 'Intergenic', 'CNV', 'Risk_freq', 'from_article', 'Accession', 'P_mlog', 'OR_Beta']
        df['tail'] = df['tail'].apply(transform_trait_id)
        df = df[['head', 'tail', 'type', 'risk allele', 'from_article']]
        df['from_article'] = df['from_article'].apply(lambda x: f'pmid{x}' if str(x).isnumeric() else x)
        df['source'] = 'gwas'
        df['label'] = 'variant_disease'
        self.data['edges']['variant_disease'] += list(df.drop_duplicates().dropna(subset=['head', 'tail']).to_dict(orient='index').values())
        return self