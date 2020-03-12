import os
import sys
import numpy as np
import networkx as nx
import scipy.io

import xml.etree.ElementTree as ET


def loadeg2Graph(filename):
    """Read .eg2 files and return the graph G,
    the number of nodes n and the number of edges m
    """
    with open(filename, "rb") as f_in:
        f_in.readline()
        G = nx.read_weighted_edgelist(f_in, nodetype=int)
    return G


def main():
    graphDirectory = sys.argv[1]
    inputDirectory = sys.argv[2]
    vectorDirectory = sys.argv[3]
    outputFilename = sys.argv[4]

    colors = np.array(['lime', 'r', 'b'])
    f_out = open(outputFilename, "w")
    for graphRoot, _, graphFiles in os.walk(graphDirectory):
        for graphFilename in graphFiles:
            dataset, ext = os.path.splitext(graphFilename)
            if ext != ".eg2":
                continue
            try:
                vectors = scipy.io.loadmat(os.path.join(vectorDirectory, "{}.mat".format(dataset)))['vec']
            except FileNotFoundError:
                print(dataset)
                continue
            graph = loadeg2Graph(os.path.join(graphRoot, graphFilename))

            for resultRoot, _, resultFiles in os.walk(inputDirectory):
                for resultFilename in resultFiles:
                    if not resultFilename.startswith(dataset):
                        continue
                    if not resultFilename.endswith("{}_graph.gexf".format(dataset)):
                        continue
                    G = nx.read_gexf(os.path.join(inputDirectory, resultFilename), int)
                    n = nx.number_of_nodes(G)
                    L = [n for n in G.nodes if 'c0' in G.nodes[n].keys()]
                    R = [n for n in G.nodes if 'c1' in G.nodes[n].keys()]
                    Lmask = np.zeros(n, dtype=bool)
                    Lmask[L] = True
                    Rmask = np.zeros(n, dtype=bool)
                    Rmask[R] = True
                    Cmask = Lmask & Rmask
                    unassigned = n - (Lmask.sum() + Rmask.sum() - Cmask.sum())
                    print("{}: Nodes={}. Unassigned={}".format(dataset, n, unassigned))

    f_out.close()


if __name__ == "__main__":
    main()

