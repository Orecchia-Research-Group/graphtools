from collections import Counter
import numpy as np
import networkx as nx
import pandas as pd

from itertools import product
np.random.seed(0)

def create_graph(nodes, cluster_prob, p, q, r, e):
    """Create a SBM graph of size `nodes` with probability of
    belonging to L, R, C according to `cluster_probs` and
    edge probability between them proportional to
    `p` for intracluster,
    `q` for cluster to overlap
    `r` for inside the overlap
    `e` for intercluster

    Each number is multiplied by log(n)/n to become a probability
    (and keep the graph connected)
    """
    clusters = np.random.choice([0, 1, 2], size=nodes, p=cluster_prob)
    cluster_sizes = Counter(clusters)
    cluster_sizes = [cluster_sizes[i] for i in sorted(cluster_sizes.keys())]

    gmean = np.sqrt(np.outer(cluster_sizes, cluster_sizes))

    probability_matrix = np.array([[p, e, q], [e, p, q], [q, q, r]]) * np.log(gmean) / gmean
    print(probability_matrix)
    G = nx.stochastic_block_model(cluster_sizes, probability_matrix)
    return G


nodes = 10000
e = 0.005
for cluster_probs in np.array([(0.49, 0.49, 0.02), (0.74, 0.24, 0.02)]):
    for p, q, r in [(1, 1, 0.2), (1, 1, 0.5), (1, 1, 1), (1, 1, 2), (1, 1, 3), (3, 2, 1), (3, 2, 2), (3, 2, 4), (3, 2, 6), (3, 4, 1), (3, 4, 3), (3, 4, 5), (3, 4, 6), (5, 3, 0.2), (10, 4, 0.4)]:
        while True:
            G = create_graph(nodes, cluster_probs, p, q, r, e)
            n = nx.number_of_nodes(G)
            m = nx.number_of_edges(G)
            conn = nx.is_connected(G)
            print(p, q, r, m, conn)
            if not conn:
                continue
            prefix = "synthetic_graph_{:.0f}_{:.0f}_{:.0f}_{:.0f}_{:.0f}".format(p*10, q*10, r*10, *cluster_probs * 100)
            with open("synthetic/{}.eg2".format(prefix), "w") as f_out:
                print(n, n, 2*m, file=f_out)
                for v in G.nodes:
                    for _, u in G.edges(v, data=False):
                        print(v, u, 1, file=f_out)
            edgelistFile = "synthetic/{}.edgelist".format(prefix)
            nx.write_edgelist(G, edgelistFile, data=False, delimiter='\t')
            break

