import gzip
import json
from collections import defaultdict, OrderedDict
import networkx as nx
import pandas as pd


categories = OrderedDict([
    ('theory', ['stoc', 'soda', 'focs', 'icalp']),
    ('dmw', ['kdd', 'icdm', 'wsdm', 'www']),
    ('ivc', ['cvpr', 'iccv', 'eccv']),
    ('ml', ['nips', 'colt', 'icml']),
    ('idm', ['sigmod', 'pods', 'vldb']),
    ('net', ['sigcomm', 'infocom', 'mobicom', 'imc']),
])

with gzip.GzipFile('dblp.json.gz', 'r') as f:
    papers = json.loads(f.read().decode())

edges = defaultdict(float)
author_areas = defaultdict(lambda: defaultdict(int))
for paper_title, paper_authors, paper_year in papers:
    for category, conferences in categories.items():
        for conference in conferences:
            if conference in paper_title.split('/'):
                for author in paper_authors:
                    author_areas[author][category] += 1
                k = len(paper_authors)
                for u in paper_authors:
                    for v in paper_authors:
                        edges[(u, v)] += 1/k


edgelist = [(u, v, int(round(w * 100))) for (u, v), w in edges.items()]
G = nx.Graph()
G.add_weighted_edges_from(edgelist)

extractedDblp = sorted(max(nx.connected_components(G), key=len))
print(type(extractedDblp))
G = G.subgraph(extractedDblp)
weights = nx.get_edge_attributes(G, 'weight')
n = len(G)
m = 2 * len(G.edges) - n
nodes = {v: i for i, v in enumerate(extractedDblp)}
with open('dcExtractedDblp.eg2', "w") as f_out:
    print(n, n, m, file=f_out)
    for v in G:
        for _, u, w in G.edges(v, data='weight'):
            print(nodes[v], nodes[u], w, file=f_out)
edgelistFile = "dcExtractedDblp.edgelist"
nx.write_edgelist(G, edgelistFile, data='weight', delimiter='\t')
areas_df = pd.DataFrame(author_areas, dtype='Int64')[extractedDblp].fillna(0).T
# print(areas_df)
areas_df.to_csv('dcExtractedDblp_authorAreas.txt', sep=' ')
