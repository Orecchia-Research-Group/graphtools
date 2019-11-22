"""Generate a planted graph


"""

import sys
import argparse
import numpy as np
import networkx as nx


def parse_args():
    """Parse arguments."""
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("nodes", type=int, help="Number of vertices in the graph.")
    parser.add_argument("--p", type=float, default=None, help="Intracluster probability of being connected.")
    parser.add_argument("--q", type=float, default=None, help="Intercluster probability of being connected.")
    parser.add_argument("--clusters", type=int, default=2, )
    parser.add_argument("--output_file", default=sys.stdout,
                        help="File to store the graph in .eg2 format. Default is standard output.")
    parser.add_argument("--seed", "-s", action="store", dest="seed", type=int, default=None, help="Randomness seed.")
    args = parser.parse_args()
    if args.p is None:
        args.p = 10 / (args.nodes / args.clusters)
    if args.q is None:
        args.q = 5 / (args.nodes / args.clusters)
    return args


def main():
    args = parse_args()
    nodes = args.nodes
    p = args.p
    q = args.q
    clusters = args.clusters
    seed = args.seed
    output_file = args.output_file

    print(nodes, p, q, clusters, output_file)

    p_array = q * np.ones([clusters, clusters])
    np.fill_diagonal(p_array, p)
    G = nx.stochastic_block_model(sizes=[nodes//clusters] * clusters, p=p_array, seed=seed)
    if not nx.is_connected(G):
        print("WARNING! Not connected")
    if isinstance(output_file, str):
        out_f = open(output_file, "w")
    else:
        out_f = output_file
    print(nodes, nodes, 2*nx.number_of_edges(G), file=out_f)
    for u, v, w in G.edges(data="weight", default=1):
        print(u, v, w, file=out_f)
        print(v, u, w, file=out_f)
    out_f.close()


if __name__ == "__main__":
    main()
