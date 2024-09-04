from neo4j import GraphDatabase
import pandas as pd
import csv

cypher = """MATCH (v1:Vocabulary)-[:MappedTo]-(:GenomicMention)-[:Mentions]-(s:Sentence)-[:Mentions]-(:GenomicMention)-[:MappedTo]-(v2:Vocabulary)
WHERE s.informative = "Informative" AND apoc.node.degree(s)<=6 AND apoc.node.degree(s)>=3 AND v1.id="{}" AND v2.id="{}"
RETURN v1.id, v2.id, s.id, s.text LIMIT 5
"""
cypher = """MATCH (v1:Vocabulary)-[r:ChemicalAffectsGeneAssociation|ChemicalOrDrugOrTreatmentToDiseaseOrPhenotypicFeatureAssociation|ChemicalToChemicalAssociation|DiseaseToPhenotypicFeatureAssociation|ExposureEventToOutcomeAssociation|GeneToDiseaseAssociation|GeneToExpressionSiteAssociation|GeneToGeneAssociation|GeneToGoTermAssociation|GeneToPathwayAssociation|VariantToDiseaseAssociation|VariantToGeneAssociation|HierarchicalStructure]-(v2:Vocabulary)
WHERE v1.id="{}" AND v2.id="{}"
RETURN v1.id, v2.id, type(r)
"""
cypher = """MATCH (s:Sentence) WHERE s.informative = "Informative" AND apoc.node.degree(s)<=6 AND apoc.node.degree(s)>=3
WITH s
MATCH (v1:ChemicalEntity|DiseaseOrPhenotypicFeature|Gene|MeshTerm)-[:MappedTo]-(:GenomicMention)-[:Mentions]-(s:Sentence)-[:Mentions]-(:GenomicMention)-[:MappedTo]-(v2:ChemicalEntity|DiseaseOrPhenotypicFeature|Gene|MeshTerm)
RETURN v1.id, v2.id, s.id, s.text
"""

def match_sents(tx):
    result = tx.run(cypher)
    return list(result)

url = "bolt://141.213.137.207:7687"

with open('/nfs/turbo/umms-drjieliu/proj/medlineKG/data/glkb_processed_data/vocab_sent_subgraph.csv', 'w', newline='') as f:
# with open('/nfs/turbo/umms-drjieliu/proj/medlineKG/data/glkb_processed_data/vocab_curated_subgraph.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows([["head","tail","sent_id","sent_text"]])
    # writer.writerows([["head","tail","rel_type"]])

    # vocabs = pd.read_pickle('/nfs/turbo/umms-drjieliu/proj/medlineKG/data/glkb_processed_data/openai_embeddings/vocab_definitions/all_vocabs_resolved.pk')
    # vocabs = vocabs[~vocabs['description'].isna()]
    # ids = vocabs['id'].to_list()
    with GraphDatabase.driver(url, auth=("neo4j", "password"), max_connection_lifetime=1000) as driver:
        with driver.session() as session:
            results = session.execute_read(match_sents)
            writer.writerows(results)
    #         for i in range(len(ids)-1):
    #             for j in range(i+1, len(ids)):
    #                 results = session.run(cypher.format(ids[i], ids[j])).values().copy()
    #                 writer.writerows(results)