import sys
import os

from collections import defaultdict

import numpy as np
import scipy.io
import scipy.sparse
import matplotlib.pyplot as plt
import networkx as nx


def read_ptn(filename):
    """Read a ptn and returns the partitions"""
    partitions = defaultdict(list)
    i = 0
    with open(filename, "r") as f_in:
        line = f_in.readline()
        if not line.startswith("cutedges"):
            for p in line.split():
                partitions[int(p)].append(i)
            i += 1
        for line in f_in:
            for p in line.split():
                partitions[int(p)].append(i)
            i += 1
    return partitions


def read_nodes(filename):
    """Read the nodes returned from Kernighan-Lin
    The nodes start from 1 and need to be reverted to start from 0"""
    nodes = []
    with open(filename, "r") as f_in:
        nodes = [int(line) - 1 for line in f_in]
    return nodes


def loadeg2Graph(filename):
    """Read .eg2 files and return the graph G,
    the number of nodes n and the number of edges m
    """
    with open(filename, "rb") as f_in:
        f_in.readline()
        G = nx.read_weighted_edgelist(f_in, nodetype=int)
    return G


def main():
    """Produce results for every graph in that directory"""
    graphDirectory = sys.argv[1]
    inputDirectory = sys.argv[2]
    outputFilename = sys.argv[3]
    colors = np.array(['lime', 'r', 'b'])
    f_out = open(outputFilename, "w")
    for graphRoot, _, graphFiles in os.walk(graphDirectory):
        for graphFilename in graphFiles:
            dataset, ext = os.path.splitext(graphFilename)
            if ext != ".eg2":
                continue
            try:
                vectors = scipy.io.loadmat(os.path.join(inputDirectory, "{}.mat".format(dataset)))['vec']
            except FileNotFoundError:
                continue
            graph = loadeg2Graph(os.path.join(graphRoot, graphFilename))

            kernighanLinPartitions = {}
            kernighanLinNodes = []
            overlap = {}
            for resultRoot, _, resultFiles in os.walk(inputDirectory):
                for resultFilename in resultFiles:
                    if not resultFilename.startswith(dataset):
                        continue
                    # Process Kernighan-Lin results
                    if resultFilename.endswith("_KernighanLin.ptn"):
                        kernighanLinPartitions = read_ptn(os.path.join(resultRoot, resultFilename))
                    elif resultFilename.endswith("_KernighanLin_nodes.txt"):
                        kernighanLinNodes = read_nodes(os.path.join(resultRoot, resultFilename))
                    elif resultFilename.endswith(".ptn") and len(resultFilename) <= len(dataset) + 9:
                        lamda = int(resultFilename.split('_')[-1].split('.')[0])
                        overlap[lamda] = read_ptn(os.path.join(resultRoot, resultFilename))
            n = nx.number_of_nodes(graph)
            m = nx.number_of_edges(graph)
            L = nx.laplacian_matrix(graph)
            overlap_perc = {}
            expansion = {}
            for lamda in sorted(overlap.keys()):
                partitions = overlap[lamda]
                L = partitions[0]
                R = partitions[1]
                Lmask = np.zeros(n, dtype=bool)
                Lmask[L] = True
                Rmask = np.zeros(n, dtype=bool)
                Rmask[R] = True
                Cmask = Lmask & Rmask
                Lonlymask = Lmask & ~Cmask
                Ronlymask = Rmask & ~Cmask
                Lonly = np.arange(n)[Lonlymask]
                Ronly = np.arange(n)[Ronlymask]
                Lsize = Lonlymask.sum()
                Rsize = Ronlymask.sum()
                colorIndex = np.zeros(n, dtype=int)
                colorIndex[Lonly] = 1
                colorIndex[Ronly] = 2
                c = colors[colorIndex]
                edgescut = len(list(nx.edge_boundary(graph, Lonly, Ronly)))
                edgescut_perc = edgescut / m
                overlapNodes = Cmask.sum()
                overlap_perc[lamda] = overlapNodes / n * 100
                expansion[lamda] = edgescut / min(Lsize, Rsize)

                plt.figure()
                ax = plt.gca()
                plt.scatter(vectors[:, 1], vectors[:, 2], c=c, s=0.25)
                plt.axis('off')
                plt.title("Graph: {}. Nodes: {}. Edges: {}. Lambda: {:d}".format(dataset, n, m, lamda))
                diff = max(vectors[:, 1]) - min(vectors[:, 1])
                ax.text(max(vectors[:, 1]) + 0.025 * diff, max(vectors[:, 2]), "Nodes: {:10d}\nEdges: {:10d}\nL: {:18d}\nR: {:18d}\nCut edges: {:d}\nOverlap: {:6d}\nExpansion: {:.3f}".format(n, m, len(L), len(R), edgescut, overlapNodes, expansion[lamda]), va='top', ha='left', fontsize=7, alpha=0.8, backgroundcolor='#d3d3d3')
                plt.savefig(os.path.join(inputDirectory, "{}_{:02d}.png".format(dataset, lamda)), bbox_inches="tight", dpi=500)
                print("Dataset: {}. Lamda: {:d}. edges: {}. overlap: {}. Expansion: {:.6f}.".format(dataset, lamda, edgescut, overlapNodes, expansion[lamda]))
                plt.close()

            if len(kernighanLinNodes):
                L = kernighanLinPartitions[0]
                R = kernighanLinPartitions[1]
                Lsize = len(L)
                Rsize = len(R)
                Lmask = np.zeros(n, dtype=bool)
                Lmask[L] = True
                Rmask = np.zeros(n, dtype=bool)
                Rmask[R] = True
                Cmask = np.zeros(n, dtype=bool)
                cutedges = len(list(nx.edge_boundary(graph, L, R)))
                KLexpansion = [cutedges / min(Lsize, Rsize)]
                KLoverlap_perc = [0]
                for i, node in enumerate(kernighanLinNodes):
                    if Lmask[node]:
                        other = np.arange(n)[Rmask & ~Cmask]
                        cutedges -= len(list(nx.edge_boundary(graph, [node], other)))
                        Rmask[node] = True
                        Rsize += 1
                    elif Rmask[node]:
                        other = np.arange(n)[Lmask & ~Cmask]
                        cutedges -= len(list(nx.edge_boundary(graph, [node], other)))
                        Lmask[node] = True
                        Lsize += 1
                    Cmask[node] = True
                    KLexpansion.append(cutedges / min(Lsize, Rsize))
                    KLoverlap_perc.append(KLoverlap_perc[-1] + 100 / n)

                lamdas = sorted(overlap.keys())
                plt.figure()
                plt.fill_between(KLoverlap_perc, KLexpansion, color='lime', alpha=0.6, label="Kernighan-Lin")
                plt.fill_between([overlap_perc[lamda] for lamda in lamdas], [expansion[lamda] for lamda in lamdas], color='blue', alpha=0.6, label="ORC-SDP")
                plt.axis(xmin=0, ymin=0)
                plt.xlabel("Overlap (%)")
                plt.ylabel("Expansion")
                plt.title("Expansion vs overlap percentage for dataset {}".format(dataset))
                plt.legend(loc="upper right")
                plt.savefig(os.path.join(inputDirectory, "over_exp_{}.png".format(dataset)), bbox_inches="tight", dpi=500)
                plt.close()

                plt.figure()
                plt.plot(lamdas, [overlap_perc[lamda] for lamda in lamdas])
                plt.axis(xmin=1, ymin=0)
                plt.xlabel(r"$\lambda$")
                plt.ylabel("Overlap (%)")
                plt.title(r"Overlap percentage vs $\lambda$ for dataset {}".format(dataset))
                plt.savefig(os.path.join(inputDirectory, "lamda_over_{}".format(dataset)), bbox_inches="tight", dpi=500)
                plt.close()


    f_out.close()


if __name__ == "__main__":
    main()
