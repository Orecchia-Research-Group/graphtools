import os
import sys
import argparse
import re

from collections import defaultdict, OrderedDict

import numpy as np
import scipy.io
import scipy.sparse
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
from sklearn.decomposition import PCA


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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
    return list(partitions.values())


def read_nodes(filename):
    """Read the nodes returned from Kernighan-Lin
    The nodes start from 1 and need to be reverted to start from 0"""
    nodes = []
    with open(filename, "r") as f_in:
        nodes = np.array([int(line) - 1 for line in f_in])
    return nodes


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


def best_partition(G, partitions):
    part = min(partitions, key=lambda p: sum(w for _, _, w in nx.edge_boundary(G, p, data="weight")))
    return part


def find_nodes(G, partition, lamda=0.5):
    n = len(G)
    partition = set(partition)
    L = partition
    degrees = np.zeros(n)
    Lsize = len(partition)
    Rsize = n - len(partition)
    Lmask = np.zeros(n, dtype=bool)
    Lmask[list(partition)] = True
    Rmask = ~Lmask
    R = np.arange(n)[Rmask]
    Cmask = np.zeros(n, dtype=bool)
    for u, v, w in nx.edge_boundary(G, partition, data="weight", default=1):
        degrees[u] += w
        degrees[v] += w
    vol = sum(d for _, d in G.degree())
    Lvol = sum(d for _, d in G.degree(L))
    Rvol = sum(d for _, d in G.degree(R))
    edgesCut = sum(degrees) / 2
    conductance = [edgesCut / min(Lvol, Rvol)]
    overlap_perc = [0]
    middle = 0
    print(edgesCut)
    while True:
        # crit = np.array([(edgesCut - degrees[v]) / min(Lvol + d, Rvol + d) for v, d in G.degree])
        idx = min(G.nodes, key=lambda v: edgesCut - degrees[v] / min(Lvol + G.degree(v), Rvol + G.degree(v)))
        maxd = degrees[idx]
        # print("Min: {:.6f} maxd: {:3.0f} degrees.max(): {:3.0f} idx: {} argmax: {}".format(crit.min(), maxd, degrees.max(), idx, degrees.argmax()))
        if maxd.max() == 0:
            break
        degrees[idx] = 0
        if Lmask[idx]:
            Rmask[idx] = True
            Cmask[idx] = True
            Rvol += G.degree(idx)
            other = Rmask & ~Cmask
            Rsize += 1
            for _, u, w in G.edges(idx, data="weight", default=1):
                if not (Rmask[u] and Lmask[u]):
                    middle += w * lamda
                if other[u]:
                    edgesCut -= w
                if Rmask[u]:
                    degrees[u] -= w
        elif Rmask[idx]:
            Lmask[idx] = True
            Cmask[idx] = True
            Lvol += G.degree(idx)
            other = Lmask & ~Cmask
            Lsize += 1
            for _, u, w in G.edges(idx, data="weight", default=1):
                if not (Rmask[u] and Lmask[u]):
                    middle += w * lamda
                if other[u]:
                    edgesCut -= w
                if Lmask[u]:
                    degrees[u] -= w
        conductance.append((edgesCut)/ min(Lvol, Rvol))
        overlap_perc.append(overlap_perc[-1] + G.degree[idx] / min(Lvol, Rvol) * 100)
    return overlap_perc, conductance


def phi(G, L, R=None, lamda=1):
    L = set(L)
    nom = 0
    Lvol = 0
    Rvol = 0
    if R is None:
        for u, v, w in G.edges(data="weight", default=1):
            Lu = u in L
            Lv = v in L
            Lvol += w * (Lu + Lv)
            if Lu ^ Lv:
                nom += w
        return nom / Lvol
    else:
        R = set(R)
        for u, v, w in G.edges(data="weight", default=1):
            Lu = u in L
            Lv = v in L
            Ru = u in R
            Rv = v in R
            if (Lu ^ Lv) and (Ru ^ Rv):
                nom += w
            Lvol += w * (Lu + Lv)
            Rvol += w * (Ru + Rv)
            if (Lu and Ru):
                nom += lamda * w
            if (Lv and Rv):
                nom += lamda * w
        return nom / min(Lvol, Rvol)


def main():
    args = parse_args()
    global verbose
    graphDirectory = args.graphDirectory
    inputDirectory = args.inputDirectory
    embedding = args.embedding
    extension = args.extension
    verbose = args.verbose

    for graphRoot, _, graphFiles in os.walk(graphDirectory):
        for graphFilename in graphFiles:
            dataset, ext = os.path.splitext(graphFilename)
            if ext != ".eg2":
                continue
            if dataset not in ["orkut","livejournal"]: continue
            # if dataset in ["ca-grQC"]: continue      # TODO: remove
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

            filenameExpr = re.compile("^{}{}_([0-9]+).ptn$".format(dataset, extension))
            nodes = {}
            partitions = {}

            if vectors.shape[1] > 2:
                pca = PCA(n_components=2)
                vectors = pca.fit_transform(vectors)
            graph = loadeg2Graph(os.path.join(graphRoot, graphFilename))
            n = len(graph)
            m = len(graph.edges)
            overlap = {}
            for resultRoot, _, resultFiles in os.walk(inputDirectory):
                for resultFilename in resultFiles:
                    if not resultFilename.startswith(dataset):
                        continue
                    if verbose:
                        print("Filename: {}".format(resultFilename))
                    match = filenameExpr.match(resultFilename)

                    # Process Sweep Cut results
                    if resultFilename.endswith("_SweepCut.ptn"):
                        partitions["Sweep Cut"] = read_ptn(os.path.join(resultRoot, resultFilename))
                    elif resultFilename.endswith("_SweepCut_nodes.txt"):
                        nodes["Sweep Cut"] = read_nodes(os.path.join(resultRoot, resultFilename))

                    # Process Sweep Cut results
                    if resultFilename.endswith("_KernighanLin.ptn"):
                        partitions["Kernighan-Lin"] = read_ptn(os.path.join(resultRoot, resultFilename))
                    elif resultFilename.endswith("_KernighanLin_nodes.txt"):
                        nodes["Kernighan-Lin"] = read_nodes(os.path.join(resultRoot, resultFilename))

                    # Process metis partition
                    elif resultFilename == "{}_metis.nodes".format(dataset):
                        nodes["METIS"] = read_nodes(os.path.join(resultRoot, resultFilename)) + 1
                    elif resultFilename == "{}_metis_1000.ptn".format(dataset):
                        partitions["METIS"] = read_ptn(os.path.join(resultRoot, resultFilename))

                    # Result of cutfind
                    elif match is not None:
                        lamda = 100 / int(match.group(1))
                        overlap[lamda] = read_ptn(os.path.join(resultRoot, resultFilename))

                    # Process bigclam
                    elif resultFilename.endswith("{}_graph.gexf".format(dataset)):
                        bigclamGraph = nx.read_gexf(os.path.join(resultRoot, resultFilename), node_type=int)
                        partitions["bigclam"] = defaultdict(list)
                        for node in bigclamGraph.nodes:
                            for cat in bigclamGraph.nodes[node].keys():
                                if not cat.startswith('c'):
                                    continue
                                partitions["bigclam"][int(cat[1:])].append(int(bigclamGraph.nodes[node]['label']))
                        partitions["bigclam"] = list(partitions["bigclam"].values())

                    # Process LEMON
                    # elif resultFilename == f"{dataset}_lemon.txt":
                    #     with open(os.path.join(resultRoot, resultFilename), "r") as lemon_in:
                    #         partitions["LEMON"] = [[int(i) for i in lemon_in]]

            # Prepare each dataset
            colors = np.array(['lime', 'red', 'blue', 'green', 'black'])
            vol = sum(d for  _, d in graph.degree())
            partition = {}
            overlap_perc, conductance = OrderedDict(), {}
            overlap_perc["ORC-SDP"] = []
            conductance["ORC-SDP"] = []
            df = pd.DataFrame(index=[r"$\phi(L)$", r"$\phi(R)$", r"$\phi(L \cup C)$", r"$\phi(R \cup C)$", r"$\phi(L, R, \lambda)$", "E(L, R)", "Overlap", "Vol(L)", "Vol(R)"])
            for lamda in sorted(overlap.keys()):
                p = overlap[lamda]
                L = p[0]
                R = p[1]
                Lmask = np.zeros(n, dtype=bool)
                Lmask[L] = True
                Rmask = np.zeros(n, dtype=bool)
                Rmask[R] = True
                Cmask = Lmask & Rmask
                Lonlymask = Lmask & ~Cmask
                Ronlymask = Rmask & ~Cmask
                Lonly = np.arange(n)[Lonlymask]
                Ronly = np.arange(n)[Ronlymask]
                C = np.arange(n)[Cmask]
                Lsize = Lonlymask.sum()
                Rsize = Ronlymask.sum()
                colorIndex = np.zeros(n, dtype=int)
                colorIndex[Lonly] = 1
                colorIndex[Ronly] = 2
                c = colors[colorIndex]

                edgescut = sum(w for _, _, w in nx.edge_boundary(graph, Lonly, Ronly, data="weight", default=1))
                vol_L = sum(d for _, d in graph.degree(L))
                vol_R = sum(d for _, d in graph.degree(R))
                vol_C = sum(d for _, d in graph.degree(C))
                overlapNodes = Cmask.sum()
                overlap_perc["ORC-SDP"].append(vol_C / min(vol_L, vol_R) * 100)
                conductance["ORC-SDP"].append((edgescut)/ min(vol_L, vol_R))
                df[lamda] = [phi(graph, Lonly), phi(graph, Ronly), phi(graph, L), phi(graph, R), phi(graph, L, R, lamda), edgescut, overlapNodes, vol_L, vol_R]
                print(f"{lamda:.2f}:", ' '.join(f"{i:.6f}" for i in df[lamda]))

                plt.figure()
                ax = plt.gca()
                plt.scatter(vectors[:, 0], vectors[:, 1], c=c, s=1, marker='o', linewidths=0)
                plt.axis('off')
                plt.title(r"Graph: {}. Nodes: {}. Edges: {}. $\lambda$: ".format(dataset, n, m) + ("1/{:.0f}" if lamda != 1 else "{:.0f}").format(1/lamda))
                diff = max(vectors[:, 0]) - min(vectors[:, 0])
                ax.text(max(vectors[:, 0]) + 0.025 * diff, max(vectors[:, 1]),
                        r"\begin{{tabular}}{{r|l}} Nodes & {} \\ Edges & {} \\ $|L|$ & {} \\ $|R|$ & {} \\ $|C|$ & {} \\ Vol($L$) & {:.0f} \\ Vol($R$) & {:.0f} \\ Vol($C$) & {:.0f} \\ Cut edges & {:.0f} \\ Conducatance & {:.6f} \end{{tabular}}"
                        .format(n, m, len(L), len(R), overlapNodes, vol_L, vol_R, vol_C, edgescut, conductance["ORC-SDP"][-1]), va='top', ha='left', fontsize=7, alpha=0.8, backgroundcolor='#d3d3d3')
                plt.savefig(os.path.join(inputDirectory, "{}_{}_k2_{:.0f}.png".format(embedding, dataset, 100 * lamda)), bbox_inches="tight", dpi=1000)
                if verbose:
                    print("Dataset: {}. Lamda: {:.2f}. edges: {}. overlap: {}. Conductance: {:.6f}.".format(dataset, lamda, edgescut, overlapNodes, conductance["ORC-SDP"][-1]))
                plt.close()
                with open(os.path.join(inputDirectory, f"{dataset}_conductance_table.txt"), "w") as f:
                    print(df.T.to_latex(float_format="%.6f"), file=f)

            for algorithm, part in sorted(partitions.items()):
                print(algorithm, end=' ')
                partition[algorithm] = best_partition(graph, part)
                overlap_perc[algorithm], conductance[algorithm] = find_nodes(graph, partition[algorithm], lamda=min(overlap.keys()))
                print(len(partition[algorithm]), len(overlap_perc[algorithm]))

            # Plot
            for algorithm in overlap_perc.keys():
                plt.plot(overlap_perc[algorithm], conductance[algorithm], label=algorithm)
                #print(algorithm, overlap_perc[algorithm])
                # print(conductance[algorithm])
            plt.axis(xmin=0, ymin=0, xmax=2*max(overlap_perc["ORC-SDP"]))
            plt.xlabel("Overlap (\%)")
            plt.ylabel("Conductance")
            plt.title(f"Conductance vs overlap percentage for dataset {dataset} and k=2")
            plt.legend(loc="upper right")
            plt.savefig(os.path.join(inputDirectory, f"over_exp_k2_{dataset}"), bbox_inches="tight", dpi=500)
            plt.close()

            # overall[dataset] = []


if __name__ == "__main__":
    main()
