import os
import sys
import networkx as nx

for root, dirs, files in os.walk(sys.argv[1]):
    for name in files:
        dataset, ext = os.path.splitext(name)
        if ext != '.eg2':
            continue
        infile = os.path.join(root, name)
        outfile = os.path.join(root, "{}.metis".format(dataset))
        if os.path.exists(outfile):
            continue
        kahip = os.path.join(root, "{}.kahip".format(dataset))
        with open(infile, "rb") as in_f:
            print(infile)
            in_f.readline()
            G = nx.read_weighted_edgelist(in_f, nodetype=int)
            G = nx.relabel_nodes(G, {n: n+1 for n in range(nx.number_of_nodes(G))})
            with open(outfile, "w") as out_f:
                print("{} {} 011".format(nx.number_of_nodes(G), nx.number_of_edges(G)), file=out_f)
                for n in sorted(G.nodes):
                    print(G.degree(n), end=' ', file=out_f)
                    print(" ".join([" ".join([str(v), str(int(w))]) for _, v, w in G.edges([n], data="weight")]), file=out_f)
            with open(kahip, "w") as out_f:
                print("{} {} 11".format(nx.number_of_nodes(G), nx.number_of_edges(G)), file=out_f)
                for n in sorted(G.nodes):
                    print(G.degree(n), end=' ', file=out_f)
                    print(" ".join([" ".join([str(v), str(int(w))]) for _, v, w in G.edges([n], data="weight")]), file=out_f)
