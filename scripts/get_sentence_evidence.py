import pandas as pd
from neo4j import GraphDatabase
from FlagEmbedding import FlagReranker
from tqdm import tqdm  # for progress tracking
from sys import argv

def init_models():
    reranker = FlagReranker('BAAI/bge-reranker-v2-m3', use_fp16=True) # Setting use_fp16 to True speeds up computation with a slight performance degradation
    return reranker

def process_sentences_batch(sentences_batch, model):
    """Process a batch of sentences and return scores"""
    pairs = [(f"Do {e1} and {e2} **directly and clearly** relate to each other in the sentence?", s) 
             for e1, e2, s, _ in sentences_batch]
    scores = model.compute_score(pairs, normalize=True)

    return scores

def process_term_pairs(term_pairs, driver, model, batch_size=1000):
    screen_query = """MATCH (v:Vocabulary {id:$id1})-[r:Cooccur]-(v2:Vocabulary {id:$id2}) RETURN r LIMIT 1"""
    cypher_query = """
    MATCH p=(v:Vocabulary {id:$id1})<-[:MappedTo]-(e1:GenomicMention)<-[:Mentions]-(s:Sentence)
    -[:Mentions]->(e2:GenomicMention)-[:MappedTo]->(v2:Vocabulary {id:$id2}) WHERE elementid(e1) <> elementid(e2)
    RETURN e1.text, e2.text, s.text, s.id LIMIT 1000
    """
    
    for _, row in tqdm(term_pairs.iterrows(), total=len(term_pairs)):
        with driver.session() as session:
            result = session.run(screen_query, id1=row['id1'], id2=row['id2']).value()
            if len(result) == 0:
                continue
            if result[0].get('evaluated'):
                continue
            # Fetch sentences
            result = session.run(cypher_query, id1=row['id1'], id2=row['id2'])
            sentences = [(r['e1.text'], r['e2.text'], r['s.text'], r['s.id']) for r in result]
            
            if not sentences:
                continue
            
            # Process sentences in batches
            positive_sent_ids = []
            for i in range(0, len(sentences), batch_size):
                batch = sentences[i:i + batch_size]
                scores = process_sentences_batch(batch, model)
                
                # Collect sentence IDs with positive scores
                batch_positive_ids = [sid for score, (_, _, _, sid) in zip(scores, batch) 
                                    if score > 0.6]
                positive_sent_ids.extend(batch_positive_ids)
            
            if positive_sent_ids:
                # Update database with evidence
                update_query = """
                MATCH (v:Vocabulary {id:$id1})-[r:Cooccur]-(v2:Vocabulary {id:$id2})
                SET r.evidence = $evidence, r.evaluate = date()
                """
                session.run(update_query, id1=row['id1'], id2=row['id2'], 
                          evidence=positive_sent_ids)
            else:
                update_query = """
                MATCH (v:Vocabulary {id:$id1})-[r:Cooccur]-(v2:Vocabulary {id:$id2})
                SET r.evaluate = date()
                """
                session.run(update_query, id1=row['id1'], id2=row['id2'])

def main():
    i = int(argv[1])
    stop_point = int(argv[2])
    bs = 1000000
    # Read term pairs
    term_pairs = pd.read_csv('/nfs/turbo/umms-drjieliu/proj/medlineKG/results/2023-05-30-create_vocab2-hyhao/fine_grained_vocab_emb/relevant_term_pair.csv')[i*bs:i*bs + bs]
    term_pairs = term_pairs.iloc[stop_point:]

    # Initialize Neo4j connection
    driver = GraphDatabase.driver(
        "bolt://141.213.137.207:7687",
        auth=("neo4j", "password"),
        max_connection_lifetime=1000
    )
    
    # Initialize models
    # tokenizer, model = init_models()
    model = init_models()
    
    try:
        process_term_pairs(term_pairs, driver, model)
    finally:
        driver.close()

if __name__ == "__main__":
    main()
