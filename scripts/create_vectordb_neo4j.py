from langchain_community.docstore.document import Document
from langchain_community.document_loaders import TextLoader
from langchain_community.vectorstores import Neo4jVector
# from langchain_openai import OpenAIEmbeddings
from langchain_text_splitters import CharacterTextSplitter
import pandas as pd
from langchain.embeddings.sentence_transformer import SentenceTransformerEmbeddings

url = "bolt://141.213.137.207:7687"
username = "neo4j"
password = "password"

# # create the open-source embedding function
model_kwargs = {'trust_remote_code': True, 'device': 'cuda'}
embedding_function = SentenceTransformerEmbeddings(model_name='Alibaba-NLP/gte-large-en-v1.5', model_kwargs=model_kwargs, multi_process=False, show_progress=False) # device_map="auto"
print('device:', embedding_function.client.device)

neo4j_vector = Neo4jVector.from_existing_graph(
url=url,
username=username,
password=password,
embedding=embedding_function,
node_label="Sentence",
embedding_node_property="embedding",
text_node_properties=["text"],
index_name='sentence_vector'
)

