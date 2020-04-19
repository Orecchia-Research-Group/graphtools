import os
import sys
import argparse
import re

from collections import defaultdict

import numpy as np
import scipy.io
import scipy.sparse
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
from sklearn.decomposition import PCA


verbose = False


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


def combine_partitions(G, partitions):
    n = nx.number_of_nodes(G)
    graph = nx.to_scipy_sparse_matrix(G)
    new_partitions = {}
    if len(partitions) > 2:
        if verbose:
            print("Combining...")
        S = None
        S_prime = None
        shared = np.zeros(n, dtype=bool)
        for _, part in sorted(partitions.items(), key=lambda x: len(x[1]), reverse=True):
            if S is None:
                S = np.zeros(n, dtype=bool)
                S[part] = True
                continue
            if S_prime is None:
                S_prime = np.zeros(n, dtype=bool)
                S_prime[part] = True
                shared |= S & S_prime
                continue
            if len(part) == 0:
                continue
            part_mask = np.zeros(n, dtype=bool)
            part_mask[part] = True
            shared |= part_mask & S
            shared |= part_mask & S_prime
            if part_mask.sum() == 0:
                continue
            S_only = np.where(S)[0]
            S_prime_only = np.where(S_prime)[0]
            part_only = np.where(part_mask & ~shared)[0]
            edges_S = len(list(nx.edge_boundary(G, part_only, S_only)))
            edges_S_prime = len(list(nx.edge_boundary(G, part_only, S_prime_only)))
            # edges_S = graph[part_mask, S & ~shared].sum()
            # edges_S_prime = graph[part_mask, S_prime & ~shared].sum()
            if edges_S > edges_S_prime:
                S |= part_mask
            else:
                S_prime |= part_mask
        new_partitions = {0: np.where(S)[0], 1: np.where(S_prime)[0]}
    else:
        new_partitions = partitions
    if verbose:
        print(', '.join([str(len(p)) for p in new_partitions.values()]))
    return new_partitions


def loadeg2Graph(filename):
    """Read .eg2 files and return the graph G,
    the number of nodes n and the number of edges m
    """
    with open(filename, "rb") as f_in:
        f_in.readline()
        G = nx.read_weighted_edgelist(f_in, nodetype=int)
    return G


def parse_args():
    parser = argparse.ArgumentParser(prog="plot")
    parser.add_argument("graphDirectory", help="Directory to read graphs from", type=str)
    parser.add_argument("inputDirectory", help="Directory to read graph partitions from", type=str)
    parser.add_argument("--embedding", "-e", choices=["spectral", "deepwalk", "node2vec", "heat_kernel"],
                        default="spectral", help="Which embedding to use for visualization")
    parser.add_argument("--extension", type=str, help="Regex to match in result files")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbosity")
    args = parser.parse_args()
    return args


def main():
    """Produce results for every graph in that directory"""
    args = parse_args()
    global verbose
    graphDirectory = args.graphDirectory
    inputDirectory = args.inputDirectory
    embedding = args.embedding
    extension = args.extension
    verbose = args.verbose

    colors = np.array(['lime', 'r', 'b'])
    bigclamGraph = None
    for graphRoot, _, graphFiles in os.walk(graphDirectory):
        for graphFilename in graphFiles:
            dataset, ext = os.path.splitext(graphFilename)
            if ext != ".eg2":
                continue
            try:
                if embedding == "spectral":
                    vectors = scipy.io.loadmat(os.path.join(inputDirectory, "{}.mat".format(dataset)))['vec'][:, 1:]
                elif embedding == "deepwalk":
                    df = pd.read_csv(os.path.join(graphDirectory, "{}_deepwalk.txt".format(dataset)), delimiter=' ', skiprows=1, index_col=0, header=None).sort_index()
                    vectors = df.values
                elif embedding == "heat_kernel":
                    vectors = scipy.io.loadmat(os.path.join(inputDirectory, "{}_heat.mat".format(dataset)))['vec']
                elif embedding == "node2vec":
                    df = pd.read_csv(os.path.join(graphDirectory, "{}_node2vec".format(dataset)), delimiter=' ', skiprows=1, index_col=0, header=None).sort_index()
                    vectors = df.values
            except FileNotFoundError:
                print("Did not find vectors for {}".format(dataset), file=sys.stderr)
                continue

            if vectors.shape[1] > 2:
                pca = PCA(n_components=2)
                vectors = pca.fit_transform(vectors)
            graph = loadeg2Graph(os.path.join(graphRoot, graphFilename))
            filenameExpr = re.compile("^{}{}_([0-9]+).ptn$".format(dataset, extension))
            kernighanLinPartitions = {}
            kernighanLinNodes = []
            overlap = {}
            for resultRoot, _, resultFiles in os.walk(inputDirectory):
                for resultFilename in resultFiles:
                    if not resultFilename.startswith(dataset):
                        continue
                    if verbose:
                        print("Filename: {}".format(resultFilename))
                    # Process Kernighan-Lin results
                    match = filenameExpr.match(resultFilename)
                    if resultFilename.endswith("_KernighanLin.ptn"):
                        kernighanLinPartitions = read_ptn(os.path.join(resultRoot, resultFilename))
                    elif resultFilename.endswith("_KernighanLin_nodes.txt"):
                        kernighanLinNodes = read_nodes(os.path.join(resultRoot, resultFilename))
                    elif match is not None:
                        lamda = int(match.group(1))
                        overlap[lamda] = read_ptn(os.path.join(resultRoot, resultFilename))
                    elif resultFilename.endswith("{}_graph.gexf".format(dataset)):
                        bigclamGraph = nx.read_gexf(os.path.join(resultRoot, resultFilename), node_type=int)
                        bigclamPartitions = defaultdict(list)
                        for node in bigclamGraph.nodes:
                            for cat in bigclamGraph.nodes[node].keys():
                                if not cat.startswith('c'):
                                    continue
                                bigclamPartitions[int(cat[1:])].append(int(bigclamGraph.nodes[node]['label']))
                        bigclamPartitions = combine_partitions(bigclamGraph, bigclamPartitions)
            for lamda, partitions in overlap.items():
                overlap[lamda] = combine_partitions(graph, partitions)
            n = nx.number_of_nodes(graph)
            m = nx.number_of_edges(graph)
            overlap_perc = {}
            expansion = {}

            L = bigclamPartitions[0]
            R = bigclamPartitions[1]
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
            edgescut = len(list(nx.edge_boundary(bigclamGraph, Lonly, Ronly)))
            edgescut_perc = edgescut / m
            overlapNodes = Cmask.sum()
            bigclamExpansion = edgescut / min(Lsize, Rsize)

            plt.figure()
            ax = plt.gca()
            plt.scatter(vectors[:, 0], vectors[:, 1], c=c, s=1, marker='o', linewidths=0)
            plt.axis('off')
            plt.title("Graph: {}. Nodes: {}. Edges: {}.".format(dataset, n, m))
            diff = max(vectors[:, 0]) - min(vectors[:, 0])
            ax.text(max(vectors[:, 0]) + 0.025 * diff, max(vectors[:, 1]), "Nodes: {:10d}\nEdges: {:10d}\nL: {:18d}\nR: {:18d}\nCut edges: {:d}\nOverlap: {:6d}\nExpansion: {:.3f}".format(n, m, len(L), len(R), edgescut, overlapNodes, bigclamExpansion), va='top', ha='left', fontsize=7, alpha=0.8, backgroundcolor='#d3d3d3')
            plt.savefig(os.path.join(inputDirectory, "{}_{}_bigclam.png".format(embedding, dataset)), bbox_inches="tight", dpi=1000)
            if verbose:
                print("Dataset: {}. bigclam edges: {}. overlap: {}. Expansion: {:.6f}.".format(dataset, edgescut, overlapNodes, bigclamExpansion))
            plt.close()

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
                plt.scatter(vectors[:, 0], vectors[:, 1], c=c, s=1, marker='o', linewidths=0)
                plt.axis('off')
                plt.title("Graph: {}. Nodes: {}. Edges: {}. Lambda: {:d}".format(dataset, n, m, lamda))
                diff = max(vectors[:, 0]) - min(vectors[:, 0])
                ax.text(max(vectors[:, 0]) + 0.025 * diff, max(vectors[:, 1]), "Nodes: {:10d}\nEdges: {:10d}\nL: {:18d}\nR: {:18d}\nCut edges: {:d}\nOverlap: {:6d}\nExpansion: {:.3f}".format(n, m, len(L), len(R), edgescut, overlapNodes, expansion[lamda]), va='top', ha='left', fontsize=7, alpha=0.8, backgroundcolor='#d3d3d3')
                plt.savefig(os.path.join(inputDirectory, "{}_{}_{:02d}.png".format(embedding, dataset, lamda)), bbox_inches="tight", dpi=1000)
                if verbose:
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




if __name__ == "__main__":
    main()
