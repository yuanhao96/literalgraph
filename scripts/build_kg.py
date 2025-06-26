import sys
from biocypher import BioCypher
from adapters.pubmed_adapter import PubmedAdapter
from adapters.journal_adapter import JournalAdapter
from adapters.vocab_adapter import OntologyAdapter, OMAdapter
from adapters.dbsnp_adapter import dbSNPAdapter
from adapters.reactome_adapter import ReactomeAdapter
from adapters.go_adapter import GOAdapter
from adapters.primekg_adapter import PrimeKGAdapter
from adapters.gwas_adapter import GWASAdapter

import os
from biocypher._logger import logger
from glob import glob
logger.debug(f"Loading module {__name__}.")

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

SCHEMA_CONFIG = "/nfs/turbo/umms-drjieliu/proj/medlineKG/data/graph_schema/glkb_schema_config.yaml"
BIOCYPHER_CONFIG = "/nfs/turbo/umms-drjieliu/proj/medlineKG/data/graph_schema/glkb_biocypher_config.yaml"

bc = BioCypher(
    biocypher_config_path=BIOCYPHER_CONFIG,
    schema_config_path=SCHEMA_CONFIG
)

files = [
    (PubmedAdapter, '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/pubmed_xml/'),
    (JournalAdapter, '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/journal_list/J_Medline.txt'),
    (OntologyAdapter, )
    (dbSNPAdapter, '/nfs/turbo/umms-drjieliu/proj/genomeKG/data/dbSNP/processed/dbSNP_snp.txt'),
    (ReactomeAdapter, {'data':'/nfs/turbo/umms-drjieliu/proj/medlineKG/data/reactome/ReactomePathways.txt', 'rt2gene':'/nfs/turbo/umms-drjieliu/proj/medlineKG/data/reactome/NCBI2Reactome.txt', 'rt2pub':'/nfs/turbo/umms-drjieliu/proj/medlineKG/data/reactome/ReactionPMIDS.txt', 'hier':'/nfs/turbo/umms-drjieliu/proj/medlineKG/data/reactome/ReactomePathwaysRelation.txt'}),
    (GOAdapter, ),
    (PrimeKGAdapter, '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/primekg/kg.csv'),
    (GWASAdapter, {'snp_gene':'/nfs/turbo/umms-drjieliu/proj/genomeKG/data/GWAS/processed/SNP_intra_gene.txt', 'snp_trait': '/nfs/turbo/umms-drjieliu/proj/genomeKG/data/GWAS/processed/SNP_trait.txt'}),
    (OMAdapter, '/nfs/turbo/umms-drjieliu/usr/xinyubao/umls_matching/database/mappings_without_dup.csv'),
]

logger.debug(bc.show_ontology_structure())

for info in files:
    if len(info) == 2: # load file from disk
        if isinstance(info[1], dict): # load multiple files at once
            adpt, p = info
            adapter = adpt().load_data(**p)
            try:
                bc.write_nodes(adapter.get_nodes(), batch_size=int(1e8))
            except StopIteration: # no nodes generated
                pass
            try:
                bc.write_edges(adapter.get_edges(), batch_size=int(1e8))
            except StopIteration: # no nodes generated
                pass
        else:
            if isinstance(info[1], list): # load list of files
                adpt, fs = info
            elif os.path.exists(info[1]):
                adpt, p = info
                if os.path.isdir(p):
                    fs = os.listdir(p)
                else:
                    fs = [p]
            
            logger.debug(f"Running {adpt.__name__}.")

            for f in fs:
                logger.debug(f"Processing data in {f}")
                adapter = adpt().load_data(os.path.join(p, f))
                
                try:
                    bc.write_nodes(adapter.get_nodes())
                except StopIteration: # no nodes generated
                    pass
                try:
                    bc.write_edges(adapter.get_edges())
                except StopIteration: # no nodes generated
                    pass

    else: # dont need to load from disk
        adpt = info[0]
        logger.debug(f"Running {adpt.__name__}.")

        adapter = adpt().load_data()
        try:
            bc.write_nodes(adapter.get_nodes(), batch_size=int(1e8))
        except StopIteration: # no nodes generated
            pass
        try:
            bc.write_edges(adapter.get_edges(), batch_size=int(1e8))
        except StopIteration: # no nodes generated
            pass

bc.summary()
bc.write_import_call()