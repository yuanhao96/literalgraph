import pyobo
import pyobo.api.utils
from collections import Counter
import gilda
from gilda import Term, make_grounder
from gilda.process import normalize
from tqdm.auto import tqdm
from biocypher._logger import logger
import json

logger.debug(f"Loading module {__name__}.")
logger.debug(pyobo.get_version())

SOURCES = {
    'gene': ["hgnc", "mesh"],
    'disease': ['mondo', "hp", "doid", "mesh"],
    'chemical': ["chebi", "mesh"],
    'anatomy': ['bto', 'uberon', 'efo', 'cl', "mesh"],
    'organism': ["ncbitaxon", "mesh"]
}

class Custom_Grounder:
    def __init__(self, prefixes:str = None, file:str = None, save_path:str = './'):
        self.save_path = save_path
        self.terms = []
        if prefixes:
            self.prefixes = prefixes
            self.term_data = []
            for prefix in self.prefixes:
                self._generate_terms(prefix)
            json.dump(self.term_data, open(self.save_path, 'w'))
        elif file:
            terms_data = json.load(open(file))
            self.terms = self.load_terms_from_file(terms_data)
        else:
            raise Exception('Please provide a list of prefixes or a file of gilda terms')
        self.grounder = gilda.make_grounder(self.terms)
    
    def load_terms_from_file(self, terms_data):
        for dat in terms_data:
            self.terms.append(
                gilda.term.Term(
                norm_text=dat.get('norm_text'), 
                text=dat.get('text'),
                db=dat.get('db'),
                id=dat.get('id'),
                entry_name=dat.get('entry_name'),
                status=dat.get('status'),
                source=dat.get('source')
                ))

    def _generate_terms(self, prefix:str):
        version = pyobo.api.utils.get_version(prefix)
        onto = pyobo.get_ontology(prefix)
        names = onto.get_id_name_mapping()
        synonyms = onto.get_id_synonyms_mapping()
        # names = pyobo.get_id_name_mapping(prefix)
        # synonyms = pyobo.get_id_synonyms_mapping(prefix)
        logger.debug(f"{prefix} v{version}, {len(names):,} names, {sum(len(v) for v in synonyms.values()):,} synonyms")

        
        for identifier, name in names.items():
            # Create a Gilda term for the standard label
            self.terms.append(gilda.Term(
                norm_text=normalize(name),
                text=name,
                db=prefix,
                id=identifier,
                entry_name=name,
                status="name",
                source=prefix,
            ))
            
            # Create a Gilda term for each synonym
            for synonym in synonyms.get(identifier, []):
                self.terms.append(gilda.Term(
                    norm_text=normalize(synonym),
                    text=synonym,
                    db=prefix,
                    id=identifier,
                    entry_name=name,
                    status="synonym",
                    source=prefix,
                ))

            # Create a Gilda term for the standard label
            self.term_data.append(gilda.Term(
                norm_text=normalize(name),
                text=name,
                db=prefix,
                id=identifier,
                entry_name=name,
                status="name",
                source=prefix,
            ).to_json())
            
            # Create a Gilda term for each synonym
            for synonym in synonyms.get(identifier, []):
                self.term_data.append(gilda.Term(
                    norm_text=normalize(synonym),
                    text=synonym,
                    db=prefix,
                    id=identifier,
                    entry_name=name,
                    status="synonym",
                    source=prefix,
                ).to_json())
            

    def ground(self, text=''):
        return self.grounder.ground(text)

    def summary(self):
        namespaces = {ns for term in self.grounder._iter_terms() for ns in term.get_namespaces()}
        status_counter = dict(Counter(term.status for term in self.grounder._iter_terms()))
        print(f"""
            Lookups: {len(self.grounder.entries):,}
            Terms: {sum(len(terms) for terms in self.grounder.entries.values())}
            Term Namespaces: {namespaces}
            Term Statuses: {status_counter}
            Adeft Disambiguators: {len(self.grounder.adeft_disambiguators)}
            """)
        if self.grounder.gilda_disambiguators:
            print(f"""Gilda Disambiguators: {len(self.grounder.gilda_disambiguators)}
            """)

class Gene_Grounder(Custom_Grounder):
    def __init__(self, prefixes:list=None, file:str=None, save_path:str='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/gene.json'):
        self.save_path = save_path
        self.terms = []
        if not prefixes:
            if file:
                self.load_terms_from_file(json.load(open(file)))
            else:
                self.term_data = []
                self.prefixes = SOURCES.get('gene')
                for prefix in self.prefixes:
                    self._generate_terms(prefix)
                json.dump(self.term_data, open(self.save_path, 'w'))
        self.grounder = gilda.make_grounder(self.terms)
        

class Disease_Grounder(Custom_Grounder):
    def __init__(self, prefixes:list = None, file:str=None, save_path:str='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/disease.json'):
        self.save_path = save_path
        self.terms = []
        if not prefixes:
            if file:
                self.load_terms_from_file(json.load(open(file)))
            else:
                self.term_data = []
                self.prefixes = SOURCES.get('disease')
                for prefix in self.prefixes:
                    self._generate_terms(prefix)
                json.dump(self.term_data, open(self.save_path, 'w'))
        self.grounder = gilda.make_grounder(self.terms)
        

class Chemical_Grounder(Custom_Grounder):
    def __init__(self, prefixes:list = None, file:str=None, save_path:str='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/chemical.json'):
        self.save_path = save_path
        self.terms = []
        if not prefixes:
            if file:
                self.load_terms_from_file(json.load(open(file)))
            else:
                self.term_data = []
                self.prefixes = SOURCES.get('chemical')
                for prefix in self.prefixes:
                    self._generate_terms(prefix)
                json.dump(self.term_data, open(self.save_path, 'w'))
        self.grounder = gilda.make_grounder(self.terms)
        

class Anatomy_Grounder(Custom_Grounder):
    def __init__(self, prefixes:list = None, file:str=None, save_path:str='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/anatomy.json'):
        self.save_path = save_path
        self.terms = []
        if not prefixes:
            if file:
                self.load_terms_from_file(json.load(open(file)))
            else:
                self.term_data = []
                self.prefixes = SOURCES.get('anatomy')
                for prefix in self.prefixes:
                    self._generate_terms(prefix)
                json.dump(self.term_data, open(self.save_path, 'w'))
        self.grounder = gilda.make_grounder(self.terms)

class Organism_Grounder(Custom_Grounder):
    def __init__(self, prefixes:list = None, file:str=None, save_path:str='/nfs/turbo/umms-drjieliu/proj/medlineKG/data/gilda_vocab/custom_grounders/organism.json'):
        self.save_path = save_path
        self.terms = []
        if not prefixes:
            if file:
                self.load_terms_from_file(json.load(open(file)))
            else:
                self.term_data = []
                self.prefixes = SOURCES.get('organism')
                for prefix in self.prefixes:
                    self._generate_terms(prefix)
                json.dump(self.term_data, open(self.save_path, 'w'))
        self.grounder = gilda.make_grounder(self.terms)

