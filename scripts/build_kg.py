import sys
from biocypher import BioCypher
sys.path.append('/home/shanghaohong/Lab/GLKB/')
from adapters.pubmed_adapter import PubmedAdapter
# from adapters.journal_adapter import JournalAdapter
# from adapters.bern2_adapter import BERN2Adapter
# from adapters.vocab_adapter import OntologyAdapter, OMAdapter
# from adapters.dbsnp_adapter import dbSNPAdapter
# from adapters.event_adapter import EventAdapter
# from adapters.semantic_relation_adapter import SemanticRelationshipAdapter
# from adapters.reactome_adapter import ReactomeAdapter
# from adapters.go_adapter import GOAdapter
# from adapters.primekg_adapter import PrimeKGAdapter
# from adapters.gwas_adapter import GWASAdapter
# from adapters.pmc_fig_adapter import PMCFigAdapter

import os
from biocypher._logger import logger
from glob import glob
logger.debug(f"Loading module {__name__}.")

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

SCHEMA_CONFIG = "configs/graph_schema/glkb_schema_config.yaml"
BIOCYPHER_CONFIG = "configs/graph_schema/glkb_biocypher_config.yaml"

bc = BioCypher(
    biocypher_config_path=BIOCYPHER_CONFIG,
    schema_config_path=SCHEMA_CONFIG
)

files = [
    (PubmedAdapter, 'files/pubmed_xml/'),
    # (JournalAdapter, '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/journal_list/J_Medline.txt'),
    # (BERN2Adapter, '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/bern2_anno/pubmed_upload/')
    # (OntologyAdapter, )
    # (dbSNPAdapter, '/nfs/turbo/umms-drjieliu/proj/genomeKG/data/dbSNP/processed/dbSNP_snp.txt'),
    # (ReactomeAdapter, {'data':'/nfs/turbo/umms-drjieliu/proj/medlineKG/data/reactome/ReactomePathways.txt', 'rt2gene':'/nfs/turbo/umms-drjieliu/proj/medlineKG/data/reactome/NCBI2Reactome.txt', 'rt2pub':'/nfs/turbo/umms-drjieliu/proj/medlineKG/data/reactome/ReactionPMIDS.txt', 'hier':'/nfs/turbo/umms-drjieliu/proj/medlineKG/data/reactome/ReactomePathwaysRelation.txt'}),
    # (GOAdapter, ),
    # (PrimeKGAdapter, '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/primekg/kg.csv'),
    # (GWASAdapter, {'snp_gene':'/nfs/turbo/umms-drjieliu/proj/genomeKG/data/GWAS/processed/SNP_intra_gene.txt', 'snp_trait': '/nfs/turbo/umms-drjieliu/proj/genomeKG/data/GWAS/processed/SNP_trait.txt'}),
    # (OMAdapter, '/nfs/turbo/umms-drjieliu/usr/xinyubao/umls_matching/database/mappings_without_dup.csv'),
    # (PMCFigAdapter, '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/figure_json_by_article/pmcimage_data.json'),
]


# p = '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/bern2_anno/pubmed_upload/'
# subs = [fld for fld in list(split(os.listdir(p), 5))]
# files = [
#     (BERN2Adapter, subs[4])
# ]
# p = '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/pubmed_xml/'
# subs = [fld for fld in list(split(os.listdir(p), 5))]
# files = [
#     (PubmedAdapter, subs[3][-1:])
# ]
# p = ''
# subs = [fld for fld in list(split(glob('/nfs/turbo/umms-drjieliu/proj/DEM_PM/DEM_PM_RES/*/*'), 5))]
# files = [
#     (EventAdapter, subs[3])
# ]
# p = '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/neo4j_csv/semantic_relationships/'
# subs = [fld for fld in list(split(os.listdir(p), 5))]
# files = [
#     (SemanticRelationshipAdapter, subs[4])
# ]

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
                adapter = adpt()
                # adapter.load_data(os.path.join(p, f))
                
                try:
                    bc.write_nodes(adapter.get_nodes(os.path.join(p, f)))
                except StopIteration: # no nodes generated
                    pass
                try:
                    bc.write_edges(adapter.get_edges(os.path.join(p, f)))
                except StopIteration: # no nodes generated
                    pass
                try:
                    pmids = adapter.get_pmids(os.path.join(p, f))
                    print(len(pmids))
                except:
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