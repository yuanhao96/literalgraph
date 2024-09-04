from langchain.embeddings.sentence_transformer import SentenceTransformerEmbeddings
# from langchain.vectorstores import Chroma
# from langchain.vectorstores import Qdrant
from qdrant_client import QdrantClient, models
from langchain.docstore.document import Document
import pandas as pd
import os
from glob import glob
import torch
import gc
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'

def pubmed2doc(fp):
    df = pd.read_csv(fp, sep=';', header=None, names=[':ID','pmcid','doi','pubmedid','title','abstract','pubdate:long','authors:string[]','author_affiliations:string[]','journal','source','publication_type:string[]','id','preferred_id',':LABEL'])
    text = (df['title'].fillna('') + ' ' + df['abstract'].fillna('')).apply(lambda x: x.strip())
    docs = [Document(page_content=doc, metadata={"source": pmid, 'level':'abstract', 'year':year}) for doc, pmid, year in zip(text, df['pubmedid'], df['pubdate:long'])]
    return docs

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def upload_data(client, collection_name, docs):
    client.upload_points(
        collection_name=collection_name,
        points=[
            models.PointStruct(
                id=doc.metadata.get('source'), vector=embedding_function.embed_query(doc.page_content), payload=doc.metadata
            )
            for doc in docs
        ],
        # shard_key_selector="user_1",
    )

# folders = '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/biocypher_glkb/biocypher-out/pubmed*'
# folders = '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/biocypher_glkb/biocypher-out/pubmed4-20240122134145/PubmedArticle-part083.csv' # /PubmedArticle-part083.csv
# folders = '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/biocypher_glkb/biocypher-out/pubmed3-20240122134118'
# folders = '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/biocypher_glkb/biocypher-out/pubmed2-20240122134050'
# folders = '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/biocypher_glkb/biocypher-out/pubmed0-20240122133952'
folders = '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/biocypher_glkb/biocypher-out/pubmed1-20240122134005'

# create the open-source embedding function
model_kwargs = {'trust_remote_code': True, 'device': 'cuda'}
embedding_function = SentenceTransformerEmbeddings(model_name='Alibaba-NLP/gte-large-en-v1.5', model_kwargs=model_kwargs, multi_process=False, show_progress=False) # device_map="auto"

client = QdrantClient(path="/nfs/turbo/umms-drjieliu/proj/medlineKG/data/vectorstores/qdrant_pubmed1")
collection_name = 'abstract'
if not client.collection_exists(collection_name):
    client.create_collection(
        collection_name=collection_name,
        vectors_config=models.VectorParams(size=1024, distance=models.Distance.COSINE, on_disk=True),
        optimizers_config=models.OptimizersConfigDiff(
            indexing_threshold=0,
        ),
        # shard_number=5,
    )
    # client.create_shard_key(collection_name, "k1")
    # client.create_shard_key(collection_name, "k2")
    # client.create_shard_key(collection_name, "k3")
    # client.create_shard_key(collection_name, "k4")
    # client.create_shard_key(collection_name, "k5")

# prev = 0
# pref = 'PubmedArticle-part'
# for folder in glob(folders):
#     # for i, fp in enumerate(glob(os.path.join(folder, f'{pref}*'))):
#     for i, fp in enumerate(glob(folder)):
#         if i < prev-1:
#             continue
#         gc.collect()
#         torch.cuda.empty_cache()
#         print('CUDA MEM usage:', torch.cuda.mem_get_info())
#         print('adding', fp)
#         docs = pubmed2doc(fp)
#         if len(docs)>0:
#             upload_data(client, collection_name, docs)

# client.update_collection(
#     collection_name=collection_name,
#     optimizer_config=models.OptimizersConfigDiff(indexing_threshold=20000),
# )

q = "SOX2 promotes lineage plasticity and antiandrogen resistance in TP53- and RB1-deficient prostate cancer"

hits = client.search(
    collection_name="abstract",
    query_vector=embedding_function.embed_query(q),
    query_filter=models.Filter(
        must=[models.FieldCondition(key="year", range=models.Range(gte=2000))]
    ),
    limit=3,
)
for hit in hits:
    print(hit.payload, "score:", hit.score)