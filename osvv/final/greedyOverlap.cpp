#include <cstdio>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <cassert>

#define min(a, b) (((a) > (b)) * (b) + ((a) <= (b)) * (a))


typedef std::map<int64_t, std::vector<std::pair<int64_t, int64_t> > > Graph;
typedef std::map<int64_t, std::set<int64_t> > Partitions;


Graph loadeg2graph(std::string graphFilename) {
    Graph G;
    std::ifstream graphFile(graphFilename);
    if (!graphFile.is_open()) {
        std::cerr << "Failed to open graphFilename " << graphFilename <<". Check path and permissions.\n";
        return G;
    }
    std::string line;
    int64_t u, v, w;
    std::getline(graphFile, line);

    while(std::getline(graphFile, line)) {
        std::istringstream iss{line};

        // Define variable 'data'. Use range constructor and stream iterator
        std::vector<int64_t> data{std::istream_iterator<int64_t>(iss), std::istream_iterator<int64_t>()};
        u = data[0];
        v = data[1];
        w = data[2];
        G[u].push_back(std::make_pair(v, w));
    }
    graphFile.close();
    return G;
}

Partitions readPtn(std::string ptnFilename) {
    Partitions partitions;
    if (!ptnFilename.size()) {
        return partitions;
    }
    std::ifstream ptnFile(ptnFilename);
    if (!ptnFile.is_open()) {
        std::cerr << "Failed to open ptnFilename " << ptnFilename <<". Check path and permissions.\n";
        return partitions;
    }
    std::string line;
    int64_t node = 0;
    while(std::getline(ptnFile, line)) {
        std::istringstream iss{line};

        // Define variable 'data'. Use range constructor and stream iterator
        std::vector<int64_t> data{std::istream_iterator<int64_t>(iss), std::istream_iterator<int64_t>()};
        for (int64_t community: data) {
            partitions[community].insert(node);
        }
        node++;
    }
    ptnFile.close();
    return partitions;
}

void greedyOverlap(
        Graph G,
        Partitions truePartitions,
        Partitions partitions,
        std::vector<int64_t> & nodes,
        std::vector<int64_t> & Lvol,
        std::vector<int64_t> & Rvol,
        std::vector<int64_t> & Cvol,
        std::vector<double> & precission,
        std::vector<double> &recall,
        std::vector<double> & f1score,
        std::vector<double> & edgeConductance,
        std::vector<double> & mixedConductance,
        std::vector<double> & lambda
    ) {

    std::map<int64_t, bool> Lmask;
    std::map<int64_t, bool> Rmask;
    std::map<int64_t, bool> Cmask;
    std::map<int64_t, int64_t> degree;
    std::map<int64_t, int64_t> edgesCut;
    std::map<int64_t, double> value;
    int64_t totalEdgesCut = 0;
    std::set<std::pair<double, int64_t> > queue;
    int64_t tp = 0;
    int64_t fp = 0;
    int64_t tr = 0;
    
    for (auto it_graph: G) {
        int64_t u = it_graph.first;
        Lmask[u] = false;
        Rmask[u] = false;
        Cmask[u] = false;
        degree[u] = 0;
        edgesCut[u] = 0;
    }

    for (auto part: partitions[0]) {
        Lmask[part] = true;
    }
    for (auto part: partitions[1]) {
        Rmask[part] = true;
    }

    nodes.push_back(-1);
    Lvol.push_back(0);
    Rvol.push_back(0);
    Cvol.push_back(0);
    precission.push_back(0);
    recall.push_back(0);
    f1score.push_back(0);
    edgeConductance.push_back(0);
    mixedConductance.push_back(0);
    lambda.push_back(1);

    for (auto it_graph: G) {
        int64_t u = it_graph.first;
        for (auto it_edge: it_graph.second) {
            int64_t v = it_edge.first;
            int64_t w = it_edge.second;
            degree[u] += w;
            if ((Lmask[u] ^ Lmask[v]) && (Rmask[u] ^ Rmask[v])) {
                edgesCut[u] += w;                
            }
        }
        totalEdgesCut += edgesCut[u];
        if (Lmask[u]) {
            Lvol[0] += degree[u];
        }
        if (Rmask[u]) { 
            Rvol[0] += degree[u];
            if (Lmask[u]) std::cerr << "Are you sure there should be nodes in the overlap at the beginning?\n";
        }
        if (truePartitions[2].find(u) != truePartitions[2].end()) {
            tr += degree[u];
        }
    }

    totalEdgesCut /= 2;     // Every edge was counted twice
    edgeConductance[0] = ((double) totalEdgesCut) / min(Lvol[0], Rvol[0]);
    mixedConductance[0] = edgeConductance[0];
    for (auto it_graph: G) {
        int64_t u = it_graph.first;
        value[u] = ((double) edgesCut[u]) / degree[u];
        queue.insert(std::make_pair(-value[u], -u));
    }

    while (totalEdgesCut > 0) {
        auto it = queue.begin();
        int64_t u = -it->second;
        lambda.push_back(-it->first);
        queue.erase(it);
        
        // fprintf(stderr, "node = %6ld edgesCut = %6ld Overlap membership = %d ", u, edgesCut[u], truePartitions[2].find(u) != truePartitions[2].end());
        nodes.push_back(u);
        totalEdgesCut -= edgesCut[u];

        edgesCut[u] = 0;

        for (auto it_edge: G[u]) {
            int64_t v = it_edge.first;
            int64_t w = it_edge.second;
            if ((Lmask[u] ^ Lmask[v]) && (Rmask[u] ^ Rmask[v])) {
                queue.erase(std::make_pair(-value[v], -v));
                edgesCut[v] -= w;
                value[v] = ((double) edgesCut[v]) / degree[v];
                queue.insert(std::make_pair(-value[v], -v));
            }
        }

        Cmask[u] = true;
        if (truePartitions[2].find(u) != truePartitions[2].end()) {
            tp += degree[u];
        } else {
            fp += degree[u];
        }
        Cvol.push_back(Cvol.back() + degree[u]);
        if (Lmask[u]) {
            Rmask[u] = true;
            Rvol.push_back(Rvol.back() + degree[u]);
            Lvol.push_back(Lvol.back());
        } else if (Rmask[u]) {
            Lmask[u] = true;
            Lvol.push_back(Lvol.back() + degree[u]);
            Rvol.push_back(Rvol.back());
        }
        
        precission.push_back(((double) tp) / tr);
        recall.push_back(((double) tp) / (tp + fp));
        if (precission.back() + recall.back() == 0) f1score.push_back(0);
        else f1score.push_back(((double) 2 * precission.back() * recall.back()) / (precission.back() + recall.back()));
        edgeConductance.push_back(((double) totalEdgesCut) / min(Lvol.back(), Rvol.back()));
        mixedConductance.push_back(((double) totalEdgesCut + lambda.back() * Cvol.back()) / min(Lvol.back(), Rvol.back()));
        // fprintf(stderr, "fp = %6ld tp = %6ld Precission = %.3lf Recall = %.3lf Conductance = %.6lf totalEdgesCut = %6ld Vol = %6ld lambda = %.6lf\n", fp, tp, precission.back(), recall.back(), conductance.back(), totalEdgesCut, min(Lvol.back(), Rvol.back()), lambda.back());
    }

}


int main(int argc, char *argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " graphFilename truthFilename ptnFilename outputFilename";
        exit(1);
    }

    Graph G = loadeg2graph(argv[1]);
    Partitions truePartitions = readPtn(argv[2]);
    Partitions partitions = readPtn(argv[3]);
    std::vector<int64_t> nodes;
    std::vector<int64_t> Lvol;
    std::vector<int64_t> Rvol;
    std::vector<int64_t> Cvol;
    std::vector<double> precission;
    std::vector<double> recall;
    std::vector<double> f1score;
    std::vector<double> edgeConductance;
    std::vector<double> mixedConductance;
    std::vector<double> lambda;
    greedyOverlap(G, truePartitions, partitions, nodes, Lvol, Rvol, Cvol, precission, recall, f1score, edgeConductance, mixedConductance, lambda);
    size_t k = nodes.size();
    assert(k == Lvol.size());
    assert(k == Rvol.size());
    assert(k == Cvol.size());
    assert(k == precission.size());
    assert(k == recall.size());
    assert(k == f1score.size());
    assert(k == edgeConductance.size());
    assert(k == mixedConductance.size());
    assert(k == lambda.size());
    std::ofstream csvFile(argv[4]);
    for (size_t i = 0; i < nodes.size(); i++) {
        csvFile << nodes[i] << ' ' << Lvol[i] << ' ' << Rvol[i] << ' ' << Cvol[i] << ' ';
        csvFile << precission[i] << ' ' << recall[i] << ' ' << f1score[i] << ' ' << edgeConductance[i] << ' ' << mixedConductance[i] << ' ' << lambda[i] << std::endl;
    }
    csvFile.close();
}

