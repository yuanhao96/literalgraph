import os
from glob import glob
import re

DATA_PATH = "/nfs/turbo/umms-drjieliu/proj/medlineKG/data/biocypher_glkb/biocypher-out"
IMPORT_PATH = "/opt/neo4j/neo4j-community-5.16.0/import/biocypher-out"
NEO4J_PATH = "/opt/neo4j/neo4j-community-5.16.0/bin/neo4j-admin"
PARAMS = [
        "--skip-bad-relationships=true",
        "--skip-duplicate-nodes=true",
        "--array-delimiter=\"|\"",
        "--delimiter=\";\"",
        "--quote='\"'",
        "--verbose",
        "--bad-tolerance=100000000000",
        # "--skip-bad-entries-logging=true",
        "--overwrite-destination=true",
        "--multiline-fields=true"
    ]
MASK=['semantic_rel']

nodes = []
edges = []
for d in glob(os.path.join(DATA_PATH, '*')):
    folder = os.path.split(d)[-1]
    if not any([folder.startswith(m) for m in MASK]):
        original_folder = folder.split('-')[-1]
        l = glob(os.path.join(d, 'neo4j-admin-import-call.sh'))
        if len(l)>0:
            import_sh = l[0]
            cmd = open(import_sh).readline()
            nodes += [s.replace(original_folder, folder) for s in re.findall("--nodes=[^\s]+", cmd)]
            edges += [s.replace(original_folder, folder) for s in re.findall("--relationships=[^\s]+", cmd)]

nodes = [n.replace(DATA_PATH, IMPORT_PATH) for n in nodes]
edges = [n.replace(DATA_PATH, IMPORT_PATH) for n in edges]
print(f"{NEO4J_PATH} database import full {' '.join(PARAMS)} {' '.join(nodes)} {' '.join(edges)}")