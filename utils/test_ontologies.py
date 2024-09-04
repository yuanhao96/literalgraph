import pyobo

ontos = ['bto', 'chebi', 'cl', 'efo', 'hgnc', 'hp', 'mondo', 'ncbitaxon', 'so', 'uberon']
for ont in ontos:
    print(ont)
    pyobo.get_ontology(ont, force=False, rewrite=False)