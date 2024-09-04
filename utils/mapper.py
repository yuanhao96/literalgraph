import pandas as pd

class biomart_mapper:
    def __init__(self, file:str = '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/ensembl_biomart/biomart220516.txt'):
        df = pd.read_csv(file, sep='\t')
        # entrez to hgnc
        sub_df = df[['HGNC ID', 'NCBI gene (formerly Entrezgene) ID']].dropna()
        entrez2hgnc = dict(zip(sub_df['NCBI gene (formerly Entrezgene) ID'].astype(int).astype(str), sub_df['HGNC ID'].astype(str).apply(lambda x: x.split(':')[-1])))
        # ensembl to hgnc
        sub_df = df[['HGNC ID', 'Gene stable ID']].dropna()
        ensembl2hgnc = dict(zip(sub_df['Gene stable ID'], sub_df['HGNC ID'].astype(str).apply(lambda x: x.split(':')[-1])))
        # self.mappings = {**self.mappings, **entrez2hgnc, **ensembl2hgnc}
        sub_df = df[['HGNC ID', 'Gene name', 'Gene Synonym']].dropna()
        name2hgnc = dict(zip(sub_df['Gene name'], sub_df['HGNC ID'].astype(str).apply(lambda x: x.split(':')[-1])))
        syn2hgnc = dict(zip(sub_df['Gene Synonym'], sub_df['HGNC ID'].astype(str).apply(lambda x: x.split(':')[-1])))
        m = {'entrez':entrez2hgnc, 'ensembl':ensembl2hgnc, 'short name':name2hgnc, 'symbol':syn2hgnc}
        self.mapper = m
    def get(self, id, namespace):
        return self.mapper.get(namespace).get(id)

class drugbank_mapper:
    def __init__(self, file:str = '/nfs/turbo/umms-drjieliu/proj/medlineKG/data/ontologies/chebi/drug-mappings.tsv'):
        mappings = pd.read_csv(file, sep='\t')
        mappings.columns = [t.strip() for t in mappings.columns]
        mappings = mappings[['drugbankId', 'chebi_id']]
        mappings['drugbankId'] = mappings['drugbankId'].str.strip()
        mappings['chebi_id'] = mappings['chebi_id'].str.strip()
        mappings = mappings[(mappings['chebi_id']!='null') & (mappings['drugbankId']!='null')]
        mappings = mappings.dropna()
        m = {'drugbank': dict(zip(mappings['drugbankId'], mappings['chebi_id']))}
        self.mapper = m
    def get(self, id, namespace='drugbank'):
        return self.mapper.get(namespace).get(id)