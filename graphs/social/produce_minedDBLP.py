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

edges = defaultdict(lambda: defaultdict(int))
author_areas = defaultdict(lambda: defaultdict(int))
for paper_title, paper_authors, paper_year in papers:
    for category, conferences in categories.items():
        for conference in conferences:
            if conference in paper_title.split('/'):
                for author in paper_authors:
                    author_areas[author][category] += 1
                for u in paper_authors:
                    for v in paper_authors:
                        if u != v:
                            edges[u][v] += 1


G = nx.from_dict_of_lists(edges)
extractedDblp = sorted(max(nx.connected_components(G), key=len))
print(type(extractedDblp))
G = G.subgraph(extractedDblp)
n = len(G)
m = 2 * len(G.edges)
nodes = {v: i for i, v in enumerate(extractedDblp)}
with open('extractedDblp.eg2', "w") as f_out:
    print(n, n, m, file=f_out)
    for v in G:
        for _, u in G.edges(v, data=False):
            print(nodes[v], nodes[u], 1, file=f_out)
edgelistFile = "extractedDblp.edgelist"
nx.write_edgelist(G, edgelistFile, data=False, delimiter='\t')
areas_df = pd.DataFrame(author_areas, dtype='Int64')[extractedDblp].fillna(0).T
print(areas_df)
areas_df.to_csv('extractedDblp_authorAreas.txt', sep=' ')
