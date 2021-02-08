import os
from collections import Counter
import numpy as np
import networkx as nx

np.random.seed(0)


def create_graph(nodes, cluster_prob, p_func, q_func, e_func):
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
    p = p_func(gmean)
    e = e_func(gmean)
    q = q_func(gmean)

    probability_matrix = np.array([
        [p[0, 0], e[0, 1], p[0, 2]],
        [e[1, 0], p[1, 1], p[1, 2]],
        [p[2, 0], p[2, 1], q[2, 2]]
    ])
    print(probability_matrix)
    G = nx.stochastic_block_model(cluster_sizes, probability_matrix)
    return G


if __name__ == '__main__':
    nodes = [10000]
    mult = 3

    e_func = lambda gmean: 0.005 * np.log(gmean) / (2 * gmean)
    for n in nodes:
        for name, cluster_probs in np.array([
                ('balanced_invSqrt', (0.50 - 1 / (2 * np.sqrt(n)), 0.50 - 1 / (2 * np.sqrt(n)), 1 / np.sqrt(n))),
                ('unbalanced_invSqrt', (0.75 - 1 / (2 * np.sqrt(n)), 0.25 - 1 / (2 * np.sqrt(n)), 1 / np.sqrt(n))),
                ('balanced_log', (0.5 - np.log(n) / (2 * n), 0.5 - np.log(n) / (2 * n), np.log(n) / (n))),
                ('unbalanced_log', (0.5 - np.log(n) / (2 * n), 0.5 - np.log(n) / (2 * n), np.log(n) / (n))),
        ]):
            cluster_probs /= np.linalg.norm(cluster_probs, 1)
            print(cluster_probs)
            for p_name, p_func in [
                    ('normal', lambda gmean: mult * np.log(gmean) / (2 * gmean)),
                    ('sqrt', lambda gmean: mult / np.sqrt(gmean)),
            ]:

                for q_name, q_func in [
                        ('sparse', lambda gmean: mult / (2 * gmean)),
                        ('normal', lambda gmean: mult * np.log(gmean) / (4 * gmean)),
                        ('dense', lambda gmean: mult / (5 * np.sqrt(gmean))),
                        ('constant', lambda gmean: mult * 1e-2 * np.ones(gmean.shape)),

                ]:
                    for k in range(5):
                        for i in range(3):
                            prefix = "synthetic_graph_{}_{}_{}_{:d}".format(name, p_name, q_name, k+1)
                            filename = "new_synthetic/{}.eg2".format(prefix)
                            if os.path.exists(filename):
                                continue;
                            G = create_graph(n, cluster_probs, p_func, q_func, e_func)
                            m = nx.number_of_edges(G)
                            conn = nx.is_connected(G)
                            print(prefix, conn)
                            if not conn:
                                continue
                            with open("new_synthetic/{}.eg2".format(prefix), "w") as f_out:
                                print(n, n, 2*m, file=f_out)
                                for v in G.nodes:
                                    for _, u in G.edges(v, data=False):
                                        print(v, u, 1, file=f_out)
                            edgelistFile = "new_synthetic/{}.edgelist".format(prefix)
                            nx.write_edgelist(G, edgelistFile, data=False, delimiter='\t')
                            break
