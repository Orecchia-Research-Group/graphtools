/* C MATLAB FUNCTION: loadmetisgraph

PURPOSE:    Reads in a hMETIS file as described in the manual 
            (http://glaros.dtc.umn.edu/gkhome/fetch/sw/hmetis/manual.pdf) for fmt = 10
            and return the graph and the node weights as specified in the file

USAGE:      function [G, weights] = loadhmetisgraph(graphFilename);

INPUTS:
-graphFilename (char):  Path to the graph file

Outputs:
-G (sparse matrix):         A sparse representation of the graph read
- weights (int64 vector):   Degree of each node

TODO: Be able to read all hMETIS formats
*/

#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"


void loadhMetisGraph(std::string graphFilename, size_t &n, size_t &m, std::vector<size_t> &heads, std::vector<size_t> &tails, std::vector<double> &weights, std::vector<double> &nodeWeights) {
    int64_t flag = 0;
    std::ifstream graphFile(graphFilename);
    if (!graphFile.is_open()) {
        std::cerr << "Failed to open graphFilename " << graphFilename <<". Check path and permissions.\n";
        return;
    }

    n = 0;
    m = 0;
    std::string line;
    size_t node = 0;

    // First line
    std::getline(graphFile, line);
    std::istringstream iss{line};
    std::vector<int64_t> data{std::istream_iterator<size_t>(iss), std::istream_iterator<size_t>()};
    m = data[0];
    n = data[1];
    if (data.size() > 2)
        flag = data[2];
    std::cerr << "Read first line" << std::endl;

    // Read hyperedges 
    for (size_t h = 0; h < m; h++) { 
        std::getline(graphFile, line);
        std::istringstream iss{line};
        std::vector<int64_t> data{std::istream_iterator<size_t>(iss), std::istream_iterator<size_t>()};

        // Internal edge
        heads.push_back(n + 2 * h + 1);
        tails.push_back(n + 2 * h);
        weights.push_back(1);

        // Star edges
        for (size_t j = 0; j < data.size(); j++) {
            heads.push_back(n + 2 * h);
            tails.push_back(data[j] - 1);
            weights.push_back(1);

            heads.push_back(data[j] - 1);
            tails.push_back(n + 2 * h + 1);
            weights.push_back(1);
        }
    }

    std::cerr << "Read hyperedges line" << std::endl;
    // Node weights
    for (size_t v = 0; v < n; v++) { 
        std::getline(graphFile, line);
        std::istringstream iss{line};
        std::vector<int64_t> data{std::istream_iterator<size_t>(iss), std::istream_iterator<size_t>()};

        nodeWeights.push_back(data[0]);
    }
    std::cerr << "Read node weights" << std::endl;

    // Hyperedge star weights are zero
    for (size_t h = 0; h < 2 * m; h++) {
        nodeWeights.push_back(0);
    }
    std::cerr << "Returning" << std::endl;
}

class MexFunction : public matlab::mex::Function
{
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    matlab::data::ArrayFactory factory;

public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {
        size_t n;
        size_t m;
        
        std::vector<size_t> heads;
        std::vector<size_t> tails;
        std::vector<double> weights;
        std::vector<double> nodeWeights;


        std::string graphFilename(matlab::data::CharArray(inputs[0]).toAscii());
        loadhMetisGraph(graphFilename, n, m, heads, tails, weights, nodeWeights);

        auto nodeWeightsArr = factory.createArray({1, n + 2 * m}, nodeWeights.begin(), nodeWeights.end());

        matlab::data::buffer_ptr_t<size_t> heads_p = factory.createBuffer<size_t>(heads.size());
        matlab::data::buffer_ptr_t<size_t> tails_p = factory.createBuffer<size_t>(tails.size());
        matlab::data::buffer_ptr_t<double> weights_p = factory.createBuffer<double>(weights.size());

        size_t *headsPtr = heads_p.get();
        size_t *tailsPtr = tails_p.get();
        double *weightsPtr = weights_p.get();


        std::for_each(heads.begin(), heads.end(), [&](const size_t& e) { *(headsPtr++) = e; });
        std::for_each(tails.begin(), tails.end(), [&](const size_t& e) { *(tailsPtr++) = e; });
        std::for_each(weights.begin(), weights.end(), [&](const double& e) { *(weightsPtr++) = e; });

        // Use the buffers to create the sparse array
        matlab::data::SparseArray<double> G =
            factory.createSparseArray<double>({n + 2 * m, n + 2 * m}, weights.size(), 
                std::move(weights_p), std::move(tails_p), std::move(heads_p));

        heads.clear();
        tails.clear();
        weights.clear();
        nodeWeights.clear();

        outputs[0] = G;
        outputs[1] = nodeWeightsArr;

    }
};
