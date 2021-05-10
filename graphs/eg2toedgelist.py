import os
import sys
import networkx as nx

for root, dirs, files in os.walk(sys.argv[1]):
    for name in files:
        dataset, ext = os.path.splitext(name)
        if ext != '.eg2':
            continue
        infile = os.path.join(root, name)
        outfile = os.path.join(root, "{}.edgelist".format(dataset))
        with open(infile, "rb") as in_f:
            print(infile)
            in_f.readline()
            G = nx.read_weighted_edgelist(in_f)
            nx.write_edgelist(G, outfile, data=False, delimiter='\t')
