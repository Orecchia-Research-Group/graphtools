from collections import Counter
import numpy as np
import networkx as nx
import pandas as pd

from itertools import product
np.random.seed(0)

def create_graph(nodes, cluster_prob, p, q, e):
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

    probability_matrix = np.array([[p, e, p], [e, q, q], [p, q, p + q]]) * np.log(gmean) / (2 * gmean)
    print(probability_matrix)
    G = nx.stochastic_block_model(cluster_sizes, probability_matrix)
    return G


nodes = [10000, 30000, 100000, 300000, 1000000]
e = 0.005
for n in nodes:
    for cluster_probs in np.array([(0.49, 0.49, 0.02), (0.74, 0.24, 0.02)]):
        for p, q in [(2, 2), (3, 2), (3, 4), (5, 3), (10, 4)]:
            for i in range(3):
                G = create_graph(n, cluster_probs, p, q, e)
                m = nx.number_of_edges(G)
                conn = nx.is_connected(G)
                print(p, q, n, m, conn)
                if not conn:
                    continue
                prefix = "synthetic_graph_{:d}_{:d}_{:d}_{:.0f}_{:.0f}".format(n, p, q, *cluster_probs * 100)
                with open("synthetic/{}.eg2".format(prefix), "w") as f_out:
                    print(n, n, 2*m, file=f_out)
                    for v in G.nodes:
                        for _, u in G.edges(v, data=False):
                            print(v, u, 1, file=f_out)
                edgelistFile = "synthetic/{}.edgelist".format(prefix)
                nx.write_edgelist(G, edgelistFile, data=False, delimiter='\t')
                break

