import openai
import pandas as pd
from neo4j import GraphDatabase

class NEREvaluator:
    def __init__(self, db_uri, db_user, db_password, openai_api_key, citation_threshold=10, fp_rate_threshold=0.5):
        self.db_uri = db_uri
        self.db_user = db_user
        self.db_password = db_password
        self.openai_api_key = openai_api_key
        self.citation_threshold = citation_threshold
        self.fp_rate_threshold = fp_rate_threshold
        openai.api_key = openai_api_key

    def get_all_vocab_ids(self):
        cypher = f"MATCH (v:Vocabulary) WHERE v.n_citation > {self.citation_threshold} RETURN v.id"
        with GraphDatabase.driver(self.db_uri, auth=(self.db_user, self.db_password)) as driver:
            with driver.session() as session:
                result = session.run(cypher)
                return [record["v.id"] for record in result]

    def get_sentences_for_vocab(self, vocab_id, limit=100):
        cypher = (
            "MATCH (s:Sentence)--(m:GenomicMention)--(v:Vocabulary) "
            "WHERE v.id = $vocab_id "
            "RETURN DISTINCT s.text AS sentence, m.text AS mention, v.id AS vocab_id, v.name AS vocab_name, labels(v) AS labels "
            "LIMIT $limit"
        )
        with GraphDatabase.driver(self.db_uri, auth=(self.db_user, self.db_password)) as driver:
            with driver.session() as session:
                result = session.run(cypher, {"vocab_id": vocab_id, "limit": limit})
                return pd.DataFrame([dict(record) for record in result])

    def evaluate_ner_with_llm(self, sentence, mention, vocab_name):
        # Prompt LLM to check if the mention is a correct extraction for the entity in the sentence
        prompt = (
            f"Sentence: \"{sentence}\"\n"
            f"NER Extraction: \"{mention}\"\n"
            f"Entity: \"{vocab_name}\"\n"
            "Is the NER extraction a correct mention of the entity in the sentence? "
            "Answer 'yes' or 'no' and briefly justify."
        )
        response = openai.ChatCompletion.create(
            model="gpt-4o-mini",
            messages=[{"role": "user", "content": prompt}],
            max_tokens=32,
            temperature=0
        )
        answer = response.choices[0].message['content'].strip().lower()
        return "yes" in answer

    def evaluate_vocab(self, vocab_id, limit=100):
        df = self.get_sentences_for_vocab(vocab_id, limit)
        if df.empty:
            return 0, 0, 0.0
        correct = 0
        for _, row in df.iterrows():
            if self.evaluate_ner_with_llm(row['sentence'], row['mention'], row['vocab_name']):
                correct += 1
        total = len(df)
        fp_rate = 1 - (correct / total) if total > 0 else 0
        return correct, total, fp_rate

    def run_evaluation(self):
        flagged_vocab_ids = []
        all_vocab_ids = self.get_all_vocab_ids()
        results = []
        for vocab_id in all_vocab_ids:
            correct, total, fp_rate = self.evaluate_vocab(vocab_id)
            results.append({
                "vocab_id": vocab_id,
                "correct": correct,
                "total": total,
                "fp_rate": fp_rate
            })
            if fp_rate > self.fp_rate_threshold:
                flagged_vocab_ids.append(vocab_id)
        return results, flagged_vocab_ids

# Example usage:
# evaluator = NEREvaluator(
#     db_uri="bolt://141.213.137.207:7687",
#     db_user="neo4j",
#     db_password="password",
#     openai_api_key="sk-proj-...",
#     citation_threshold=10,
#     fp_rate_threshold=0.5
# )
# results, flagged = evaluator.run_evaluation()
# print("Flagged vocabularies:", flagged)
