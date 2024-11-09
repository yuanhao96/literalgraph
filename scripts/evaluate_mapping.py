from langchain_openai import ChatOpenAI, OpenAIEmbeddings
from langchain_core.prompts import ChatPromptTemplate
from langchain_community.docstore.document import Document
from langchain_community.document_loaders import TextLoader
from langchain_community.vectorstores import Neo4jVector
from langchain.embeddings.sentence_transformer import SentenceTransformerEmbeddings
from neo4j import GraphDatabase
import pandas as pd
from difflib import SequenceMatcher
from collections import Counter

openai_key = ""
glkb_bolt_url = "bolt://141.213.137.207:7687"
llm = "gpt-4o-mini"
model_path = 'Alibaba-NLP/gte-large-en-v1.5'

model_kwargs = {'trust_remote_code': True}
embedding_function = SentenceTransformerEmbeddings(model_name=model_path, model_kwargs=model_kwargs) # device_map="auto"

llm = ChatOpenAI(
    model=llm, 
    temperature=0,
    api_key=openai_key,
)
system = """You are an expert at biomedical question answering. \

Task:Answer biomedical questions in one sentence based on the contexts.
Instructions:
Don't try to make up an answer, if you don't know just say that you don't know.
Use only the following pieces of context to answer the question.
The context are dictionaries of PubMed articles, containing their titles, abstracts, and PubMed IDs.

Note: 
If the question is not related to the contexts, just say that it is impossible to answer the question based on the contexts.
Do not include any text except the generated answer to the question.
"""
prompt = ChatPromptTemplate.from_messages(
    [
        ("system", system),
        ("human", "question: {question} context: {context}"),
    ]
)
answer_generator = prompt | llm

retrieval_query = """
                RETURN node {.pubmedid, .title, .abstract} AS text, score, {} AS metadata
                """
abstract_store = Neo4jVector.from_existing_index(
    embedding_function,
    url=glkb_bolt_url,
    username="neo4j",
    password="password",
    index_name='abstract_vector',
    retrieval_query=retrieval_query,
)

driver = GraphDatabase.driver(glkb_bolt_url, auth=("neo4j", "password"), max_connection_lifetime=1000)

def retrieve_abstract(query, k=10):
    retrieved_docs = abstract_store.similarity_search(query, k=k)
    result = []
    for res in retrieved_docs:
        abstract, title, pmid = "", "", ""
        for r in res.page_content.strip().split('\n'):
            if r.startswith('abstract'):
                abstract = ': '.join(r.split(': ')[1:])
            elif r.startswith('title'):
                title = ': '.join(r.split(': ')[1:])
            elif r.startswith(('pubmedid')):
                pmid = ': '.join(r.split(': ')[1:])
        result.append([pmid, title, abstract])
    return result

def generate_rag_answer(question, context):
    """
    Generate answer using RAG on retrieved documents

    Returns:
        state (dict): New key added to state, generation, that contains LLM generation
    """

    # RAG generation
    generation = answer_generator.invoke({"context": context, "question": question}).content
    
    return generation

def generate_description(name, k=3):
    qs = {
        "def": "what is {}?",
        "func": "what is the biological function of of {}?",
        "gene": "genes or pathways related to {}?",
        "disease": "diseases or phenotypes related to {}?",
        "chems": "drugs or small molecule chemicals related to {}?",
    }
    res = []

    for q in qs.values():
        contexts = retrieve_abstract(q.format(name), k=k)
        ans = generate_rag_answer(q.format(name), contexts)
        if ans == "It is impossible to answer the question based on the contexts.":
            continue
        res.append(ans)
    return ' '.join(res)

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

system = """You are an expert at biomedical question answering. \

Task:classify if an extracted entity refers to a standardized biomedical term from a list of contexts that contain the entity.
Instructions:
Answer in a list of 'yes' or 'no'.
The context are dictionaries of PubMed abstracts that contain the entity.

Example output: ['yes', 'no', 'yes']

Note: 
Don't try to make up an answer, if you don't know, return 'unknown'.
Do not include any text except the generated yes/no/unknown answer.
"""
prompt = ChatPromptTemplate.from_messages(
    [
        ("system", system),
        ("human", "does the entity '{normalized_name}' refer to '{label}' '{name}' in the following contexts? context: {context}"),
    ]
)
mapping_classifier = prompt | llm

def filter_score(_id):
    if _id.startswith('doid') or _id.startswith('mondo') or _id.startswith('hp'): # disease
        return 0.8
    elif _id.startswith('chebi'): # chem
        return 0.8
    elif _id.startswith('hgnc'): # gene
        return 1.0
    else:
        return 1.0

def get_entity_context(_id, normalized_name):
    with driver.session() as session:
        res = session.run("match (v:Gene|DiseaseOrPhenotypicFeature|ChemicalEntity|MeshTerm)<-[:MappedTo]-(e:GenomicMention)<--(a:Article) where v.id='{}' and e.text=\"{}\" return distinct a.abstract limit 100".format(_id, normalized_name)).value().copy()
    if len(res) > 0:
        return res
    else:
        return None

def get_label(_id):
    if _id.startswith('doid') or _id.startswith('mondo') or _id.startswith('hp'): # disease
        return 'disease'
    elif _id.startswith('chebi'): # chem
        return 'chemical'
    elif _id.startswith('hgnc'): # gene
        return 'gene'
    elif _id.startswith('mesh'): # mesh
        return 'MeSH Term'
    return ''

def true_map(pred_mapped):
    c = Counter(pred_mapped)
    if c['yes'] < 6 or c['no'] > 0:
        return False
    return True

df = pd.read_csv('/nfs/turbo/umms-drjieliu/proj/medlineKG/results/2023-05-30-create_vocab2-hyhao/fine_grained_vocab_emb/filtered_vocab.csv')
sub = df[df['score'] < df['id'].apply(filter_score)]
sub['contexts_of_normalized_name'] = sub.apply(lambda x: get_entity_context(x['id'], x['normalized_name']), axis=1)
sub['label'] = sub['id'].apply(get_label)
# sub = pd.read_csv('/nfs/turbo/umms-drjieliu/proj/medlineKG/results/2023-05-30-create_vocab2-hyhao/fine_grained_vocab_emb/unmatched_vocab.csv')
sub['contexts_of_normalized_name'] = sub['contexts_of_normalized_name'].apply(eval)
sub.loc[sub['pred_mapped'].isna(), 'pred_mapped'] = sub[sub['pred_mapped'].isna()].apply(lambda x: mapping_classifier.invoke({
    'name': x['name'],
    'label': x['label'],
    'normalized_name': x['normalized_name'],
    'context': x['contexts_of_normalized_name']
}).content, axis=1)
sub.to_csv('/nfs/turbo/umms-drjieliu/proj/medlineKG/results/2023-05-30-create_vocab2-hyhao/fine_grained_vocab_emb/unmatched_vocab.csv', index=False)
sub[sub['pred_mapped'].apply(true_map)][['name','id','normalized_name','n_citation','score']].to_csv('/nfs/turbo/umms-drjieliu/proj/medlineKG/results/2023-05-30-create_vocab2-hyhao/fine_grained_vocab_emb/filtered_unmatched_vocab.csv', index=False)
sub[~sub['pred_mapped'].apply(true_map)][['name','id','normalized_name','n_citation','score']].to_csv('/nfs/turbo/umms-drjieliu/proj/medlineKG/results/2023-05-30-create_vocab2-hyhao/fine_grained_vocab_emb/blacklisted_unmatched_vocab.csv', index=False)