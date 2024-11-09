import pandas as pd
from neo4j import GraphDatabase

glkb_bolt_url = "bolt://141.213.137.207:7687"
driver = GraphDatabase.driver(glkb_bolt_url, auth=("neo4j", "password"), max_connection_lifetime=1000)

def remove_blacklisted_mappings(_id):
    with driver.session() as session:
        res = session.run(f"""call apoc.periodic.iterate('MATCH (v:Vocabulary)<-[r:ContainTerm]-(:Article) WHERE v.id="{_id}" AND (r.source="BERN2" OR r.source="DEM") return r', "delete r", {{batchSize:1, parallel:true}})""")
        res = session.run(f"""call apoc.periodic.iterate('MATCH (v:Vocabulary)<-[r:MappedTo]-(:GenomicMention) WHERE v.id="{_id}" return r', "delete r", {{batchSize:1, parallel:true}})""")
        res = session.run(f"""call apoc.periodic.iterate('MATCH (v:Vocabulary)-[r:Associate|Bind|Comparison|Cotreatment|NegativeCorrelation|PositiveCorrelation]-(v2:Vocabulary) WHERE v.id="{_id}" return r', "delete r", {{batchSize:1, parallel:true}})""")
        # res = session.run(f'match (v:Vocabulary)<-[r:ContainTerm]-(:Article) where v.id="{_id}" and (r.source="BERN2" or r.source="DEM") delete r')
        # res = session.run(f'match (v:Vocabulary)<-[r:MappedTo]-(:GenomicMention) where v.id="{_id}" delete r')
        # res = session.run(f'match (v:Vocabulary)-[r:Associate|Bind|Comparison|Cotreatment|NegativeCorrelation|PositiveCorrelation]-(v2:Vocabulary) where v.id="{_id}" delete r')

blacklisted = pd.read_csv('/nfs/turbo/umms-drjieliu/proj/medlineKG/results/2023-05-30-create_vocab2-hyhao/fine_grained_vocab_emb/blacklisted_unmatched_vocab.csv')
blacklisted = blacklisted[blacklisted['id'].str.contains('mesh')]

print(len(blacklisted))
i = 0
for _id in blacklisted['id']:
    print(_id)
    i += 1
    print(i, len(blacklisted))
    remove_blacklisted_mappings(_id)