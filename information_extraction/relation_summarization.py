from Graph.neo4j import run_cypher
from NLP.re import RE_SCHEMA
import json
from typing import List, Dict, Tuple, Optional
from openai import OpenAI
from neo4j import GraphDatabase

RE_SCHEMA = {
    'disease': {
        'variant': ['positve correlation', 'negative correlation', 'association'],
        'gene': ['positve correlation', 'negative correlation', 'association'],
        'chemical': ['positve correlation', 'negative correlation', 'association'],
    },
    'chemical': {
        'chemical': ['positve correlation', 'negative correlation', 'association', 'bind', 'drug interaction', 'cotreatment', 'comparison', 'conversion'],
        'variant': ['positve correlation', 'negative correlation', 'association'],
    },
    'variant': {
        'variant': ['association'],
    },
    'gene': {
        'chemical': ['positive correlation', 'negative correlation', 'association', 'bind'],
    },
    'gene': {
            'gene': ['positive correlation', 'negative correlation', 'association', 'bind'],
        },
    }

class RelationSummarizer:
    def __init__(self, db_uri, db_user, db_password, openai_api_key, re_schema, model="gpt-4o-mini-2024-07-18"):
        self.db_uri = db_uri
        self.db_user = db_user
        self.db_password = db_password
        self.openai_api_key = openai_api_key
        self.re_schema = re_schema
        self.model = model
        self.client = OpenAI(api_key=openai_api_key)

    def get_contexts_for_term_pair(self, ids1: List[str], ids2: List[str], batch_size=30) -> List[Dict]:
        node_citation_cypher = "MATCH (v:Vocabulary) WHERE elementId(v) IN $ids RETURN v.n_citation"
        cypher_query = (
            "MATCH p=(v:Vocabulary)<--(s:Article)-->(v2:Vocabulary) "
            "WHERE elementId(v) IN $ids1 AND elementId(v2) IN $ids2 AND elementId(v) <> elementId(v2) "
            "RETURN s.pubmedid AS pubmedid, s.title AS title, s.abstract AS abstract "
            "ORDER BY s.n_citation DESC LIMIT $batch_size"
        )
        cypher_query2 = (
            "MATCH p=(v:Vocabulary)<--(s:Article)-->(v2:Vocabulary) "
            "WHERE elementId(v) IN $ids1 AND elementId(v2) IN $ids2 AND elementId(v) <> elementId(v2) "
            "RETURN s.pubmedid AS pubmedid, s.title AS title, s.abstract AS abstract LIMIT $batch_size"
        )
        with GraphDatabase.driver(self.db_uri, auth=(self.db_user, self.db_password)) as driver:
            with driver.session() as session:
                result = session.run(node_citation_cypher, {"ids": ids1 + ids2}).value()
                use_simple_query = any(r is not None and r > 2000 for r in result)
                cypher = cypher_query2 if use_simple_query else cypher_query
                result = session.run(cypher, {"ids1": ids1, "ids2": ids2, "batch_size": batch_size})
                sentences = []
                sids = set()
                for r in result:
                    if r['pubmedid'] not in sids:
                        sentences.append({'title': r['title'], 'abstract': r['abstract'], 'pubmedid': r['pubmedid']})
                        sids.add(r['pubmedid'])
        return sentences

    def create_summary_prompt(self, entities: Tuple[str, ...], contexts: List[Dict]) -> str:
        combined_contexts = "\n\n".join([
            f"Article PubMed ID: {ctx['pubmedid']}\nTitle: {ctx['title']}\nAbstract: {ctx['abstract']}"
            for ctx in contexts
        ])
        prompt = f"""Analyze the following scientific articles and summarize the most significant relationship between {' and '.join(entities)}. 
If the entities are not related, set the relationship to 'co-occurring', summary to '<entity1> and <entity2> co-occur in the provided articles', and pubmedids to the pubmedids of articles that contain both entities.
If the entities do not present in any provided article, set the relationship to 'none', summary to '<entity1> and <entity2> do not present in the provided articles', and pubmedids to [].

Articles:
{combined_contexts}

Schema of possible relationships between entity types:
{json.dumps(self.re_schema, indent=2)}

Provide a summary in JSON format:
{{
    "entity1": "entity_name",
    "entity2": "entity_name",
    "entity1_type": "entity_type",
    "entity2_type": "entity_type",
    "relationship": "relationship_type",
    "summary": "summarize the relationship in a short sentence, e.g., entity1 activates entity2 via phosphorylation",
    "pubmedids": [<supporting article pubmedid>, <supporting article pubmedid>, ...],
}}

Focus on the most significant and well-supported relationship between the target entities across all articles."""
        return prompt

    def summarize_relationships(self, entities: Tuple[str, ...], contexts: List[Dict]) -> Dict:
        prompt = self.create_summary_prompt(entities, contexts)
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": "You are a scientific relationship extraction assistant. Summarize the most significant relationships between entities from scientific literature accurately and precisely."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0,
                response_format={"type": "json_object"}
            )
            result = json.loads(response.choices[0].message.content)
            # Validate result
            for key in ['entity1', 'entity2', 'relationship', 'summary', 'pubmedids']:
                if key not in result:
                    raise ValueError(f"Missing key '{key}' in LLM response")
            if not isinstance(result['pubmedids'], list):
                result['pubmedids'] = []
            return result
        except Exception as e:
            print(f"Error processing contexts: {e}")
            return {"relationships": []}

    def set_relationship_in_db(self, id1, id2, rel: Dict):
        cypher = (
            "MATCH (v1)-[r:Cooccur]-(v2) "
            "WHERE elementId(v1)=$id1 AND elementId(v2)=$id2 "
            "SET r.relationship=$relationship, r.summary=$summary, r.pubmedids=$pubmedids "
            "RETURN r"
        )
        with GraphDatabase.driver(self.db_uri, auth=(self.db_user, self.db_password)) as driver:
            with driver.session() as session:
                session.run(
                    cypher,
                    {
                        "id1": id1,
                        "id2": id2,
                        "relationship": rel.get('relationship'),
                        "summary": rel.get('summary'),
                        "pubmedids": rel.get('pubmedids', [])
                    }
                )

    def process_and_set_relationships(self, id1, id2, entities, batch_size=30):
        contexts = self.get_contexts_for_term_pair([id1], [id2], batch_size)
        rel = self.summarize_relationships(entities, contexts)
        if rel.get('relationship'):
            self.set_relationship_in_db(id1, id2, rel)
        return rel

# Example usage:
# summarizer = RelationSummarizer(
#     db_uri="bolt://141.213.137.207:7687",
#     db_user="neo4j",
#     db_password="password",
#     openai_api_key="sk-proj-...",
#     re_schema=RE_SCHEMA
# )
# result = summarizer.process_and_set_relationships(id1, id2, (entity1_name, entity2_name))
# print(result)