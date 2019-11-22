"""Generate a powerlaw cluster graph


"""

import sys
import argparse
import networkx as nx


def parse_args():
    """Parse arguments."""
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("nodes", type=int, help="Number of vertices in the graph.")
    parser.add_argument("edges", type=int, help="Number of random edges to add for each new node.")
    parser.add_argument("-p", type=float, default=0, help="Probability of adding a triangle after adding a random edge.")
    parser.add_argument("--output_file", default=sys.stdout,
                        help="File to store the graph in .eg2 format. Default is standard output.")
    parser.add_argument("--seed", "-s", action="store", dest="seed", type=int, default=None, help="Randomness seed.")
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    nodes = args.nodes
    edges = args.edges
    p = args.p
    seed = args.seed
    output_file = args.output_file

    # print(nodes, radius, dim, p, seed, output_file)

    G = nx.powerlaw_cluster_graph(n=nodes, m=edges, p=p, seed=seed)
    if not nx.is_connected(G):
        print("WARNING! Not connected")
    if isinstance(output_file, str):
        out_f = open(output_file, "w")
    else:
        out_f = output_file
    # out_f.write("{} {} {}\n".format(nodes, nodes, nx.number_of_edges(G)).encode("utf-8"))
    print(nodes, nodes, 2*nx.number_of_edges(G), file=out_f)
    for u, v, w in G.edges(data="weight", default=1):
        print(u, v, w, file=out_f)
        print(v, u, w, file=out_f)
    # nx.write_weighted_edgelist(G, out_f)
    out_f.close()


if __name__ == "__main__":
    main()
