import json
from enum import Enum
from adapters import Adapter, Node, Edge
from utils import escape_text

class NERAdapter_NodeType(Enum):
    MENTION = "genomic_mention"

class NERAdapter_Mention_Field(Enum):
    TEXT = "text"
    CHAR_START = "char_start"
    CHAR_END = "char_end"
    TYPE = "type"
    SCORE = "score"
    GROUNDINGS = "groundings"

class NERAdapter_MapTo_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for mapped_to edges.
    """
    SOURCE = "source"
    SCORE = "score"

class NERAdapter_ContainMention_EdgeField(Enum):
    """
    Define possible fields the adapter can provide for contain_mention edges.
    """
    SOURCE = "source"

class NERAdapter_ContainTerm_EdgeField(Enum):
    TEXT = "text"
    CHAR_START = "char_start"
    CHAR_END = "char_end"
    TYPE = "type"
    SOURCE = "source"
    NORM_NAME = "normalized_name"
    PROB = "prob"

class NERAdapter_EdgeType(Enum):
    CONTAIN_MENTION = "contain_mention"
    MAP_TO = "mapped_to"

class NERAdapter(Adapter):
    def __init__(self, node_types=None, node_fields=None, edge_types=None, edge_fields=None):
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)
        self.data = None

    def load_data(self, file: str):
        with open(file, "r") as f:
            self.data = json.load(f)

    def get_nodes(self, file: str = None):
        if file:
            self.load_data(file)
        nodes = []
        node_id = 0
        for entity_type, entities in self.data.items():
            for ent in entities:
                node = GenomicMention(
                    id=f"mention_{node_id}",
                    fields=[
                        NERAdapter_Mention_Field.TEXT.value,
                        NERAdapter_Mention_Field.CHAR_START.value,
                        NERAdapter_Mention_Field.CHAR_END.value,
                        NERAdapter_Mention_Field.TYPE.value,
                        NERAdapter_Mention_Field.SCORE.value,
                        NERAdapter_Mention_Field.GROUNDINGS.value
                    ],
                    properties={
                        "text": ent["entity"],
                        "char_start": ent["start"],
                        "char_end": ent["end"],
                        "type": entity_type,
                        "score": ent["score"],
                        "groundings": ent.get("groundings", [])
                    }
                )
                nodes.append(node)
                node_id += 1
        return nodes

    def get_edges(self, file: str = None):
        if file:
            self.load_data(file)
        edges = []
        node_id = 0
        for entity_type, entities in self.data.items():
            for ent in entities:
                mention_id = f"mention_{node_id}"
                # ContainMention edge (from a dummy doc node "doc_0" to mention)
                contain_edge = ContainMention(
                    source="doc_0",
                    target=mention_id,
                    id=f"contain_{node_id}",
                    fields=["source", "target"],
                    properties={}
                )
                edges.append(contain_edge)
                # MapTo edges for each grounding
                for idx, grounding in enumerate(ent.get("groundings", [])):
                    curie, score = grounding
                    mapto_edge = MapTo(
                        source=mention_id,
                        target=curie,
                        id=f"mapto_{node_id}_{idx}",
                        fields=["source", "target", "score"],
                        properties={"score": score}
                    )
                    edges.append(mapto_edge)
                    contain_term_edge = ContainTerm(
                        source="doc_0",
                        target=curie,
                        id=f"containterm_{node_id}_{idx}",
                        properties={
                            "text": ent["entity"],
                            "char_start": ent["start"],
                            "char_end": ent["end"],
                            "type": entity_type,
                            "source": "NER",
                            "norm_name": curie,
                            "prob": ent.get("score", None)
                        }
                    )
                    edges.append(contain_term_edge)
                node_id += 1
        return edges


class GenomicMention(Node):
    """
    genomic mention nodes
    """
    def __init__(self, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.label = 'genomic_mention'
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in NERAdapter_Mention_Field]
        else:
            return fields

    def _generate_properties(self, properties):
        prop_dict = {}
        for field in self.fields:
            f = properties.get(field)
            if f:
                if isinstance(f, str):
                    prop_dict[field] = escape_text(f)
                else:
                    prop_dict[field] = properties.get(field)
        
        return prop_dict

class ContainTerm(Edge):
    """
    contain_term edges
    """
    def __init__(self, source:str, target:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = "contain_term"
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in NERAdapter_ContainTerm_EdgeField]
        else:
            return fields

    def _generate_properties(self, properties):
        prop_dict = {}
        for field in self.fields:
            f = properties.get(field)
            if f:
                if isinstance(f, str):
                    prop_dict[field] = escape_text(f)
                else:
                    prop_dict[field] = properties.get(field)
        
        return prop_dict

class ContainMention(Edge):
    """
    contain_mention edges
    """
    def __init__(self, source:str, target:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = "contain_mention"
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in NERAdapter_EdgeType]
        else:
            return fields

    def _generate_properties(self, properties):
        prop_dict = {}
        for field in self.fields:
            f = properties.get(field)
            if f:
                if isinstance(f, str):
                    prop_dict[field] = escape_text(f)
                else:
                    prop_dict[field] = properties.get(field)
        
        return prop_dict


class MapTo(Edge):
    """
    mapped_to edges
    """
    def __init__(self, source:str, target:str, id:str=None, fields: list = None, properties: dict = None):
        self.id = id
        self.source = source
        self.target = target
        self.label = "mapped_to"
        self.fields = self._generate_fields(fields)
        self.properties = self._generate_properties(properties)
    
    def _generate_fields(self, fields):
        if fields is None:
            return [i.value for i in NERAdapter_MapTo_EdgeField]
        else:
            return fields

    def _generate_properties(self, properties):
        prop_dict = {}
        for field in self.fields:
            f = properties.get(field)
            if f:
                if isinstance(f, str):
                    prop_dict[field] = escape_text(f)
                else:
                    prop_dict[field] = properties.get(field)
        
        return prop_dict

# Example usage
if __name__ == "__main__":
    adapter = NERAdapter()
    adapter.load_data("example_ner_output.json")
    nodes = adapter.get_nodes()
    edges = adapter.get_edges()
    print("Nodes:")
    for n in nodes:
        print(vars(n))
    print("\nEdges:")
    for e in edges:
        print(vars(e))