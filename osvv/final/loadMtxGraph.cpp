/* C MATLAB FUNCTION: loadMtxGraph

PURPOSE:    Reads in a Matrix Market file.

USAGE:      function [G, weights] = loadMtxGraph(graphFilename);

INPUTS:
- graphFilename (char):  Path to the graph file

Outputs:
- G (sparse matrix):         A sparse representation of the graph read
- weights (int64 vector):   Degree of each node

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


bool verifyBanner(std::string graphFilename, std::string line, bool &symmetric) {
    std::istringstream iss{line};
    std::string data[5];
    iss >> data[0] >> data[1] >> data[2] >> data[3] >> data[4];
    // std::vector<std::string> data{std::istream_iterator<size_t>(iss), std::istream_iterator<size_t>()};
    if (data[0] != "%%MatrixMarket") {
        std::cerr << "Matrix Market Banner not present in " << graphFilename << ".\n";
        return false;
    }

    if (data[1] != "matrix") {
        std::cerr << "You have a dense matrix in " << graphFilename << ". Are you sure this is a good idea?.\n";
        return false;
    }

    if (data[2] != "coordinate") {
        std::cerr << "Somehow you don't have coordinates in " << graphFilename << ". Not sure what to do with that.\n";
        return false;
    }

    if (data[3] != "integer") {
        std::cerr << "Flow for the cut/matching game works only for integer weights. " << graphFilename << " has " << data[3] << " type weights and is not supported.\n";
        return false;
    }

    if (data[4] == "general") {
        symmetric = false;
    } else if (data[4] == "symmetric") {
        symmetric = true;
    } else {
        std::cerr << "No idea what to do with hermitian or skew-symmetric in " << graphFilename <<".\n";
        return false;
    }
    return true;
}


void loadhMetisGraph(std::string graphFilename, size_t &n, size_t &m, std::vector<size_t> &heads, std::vector<size_t> &tails, std::vector<double> &weights, std::vector<double> &nodeWeights) {
    int64_t flag = 0;
    std::string suffix = ".mtx";
    if (!graphFilename.ends_with(suffix)) {
        throw std::runtime_error("File needs to be .mtx");
    }

    std::ifstream graphFile(graphFilename);
    if (!graphFile.is_open()) {
        std::cerr << "Failed to open graphFilename " << graphFilename <<". Check path and permissions.\n";
        return;
    }

    n = 0;
    m = 0;
    bool symmetric = false;
    std::string line;
    size_t node = 0;

    // Verify banner and detect if symmetric
    std::getline(graphFile, line);
    if (!verifyBanner(graphFilename, line, symmetric)) {
        return;
    }

    while (graphFile.peek() == '%') graphFile.ignore(2048, '\n');

    // First line
    std::getline(graphFile, line);
    std::istringstream iss{line};
    std::vector<int64_t> data{std::istream_iterator<size_t>(iss), std::istream_iterator<size_t>()};
    if (data[0] != data[1]) {
        std::cerr << "Only able to process square adjacency lists. " << graphFilename << " has sizes " << data[0] << " and " << data[1] <<".\n";
        return;
    }
    n = data[0];
    m = data[2];

    for (size_t i = 0; i < n; i++)
        nodeWeights.push_back(0);

    // Read edges
    for (size_t h = 0; h < m; h++) {
        std::getline(graphFile, line);
        std::istringstream iss{line};
        std::vector<int64_t> data{std::istream_iterator<size_t>(iss), std::istream_iterator<size_t>()};

        // Internal edge
        tails.push_back(data[0] - 1);
        heads.push_back(data[1] - 1);
        weights.push_back(data[2]);
        nodeWeights[data[0] - 1] += data[2];

        if (symmetric) {
            tails.push_back(data[1] - 1);
            heads.push_back(data[0] - 1);
            weights.push_back(data[2]);
            nodeWeights[data[1] - 1] += data[2];
        }
    }


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

        auto nodeWeightsArr = factory.createArray({1, n}, nodeWeights.begin(), nodeWeights.end());

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
                factory.createSparseArray<double>({n, n}, weights.size(),
                                                  std::move(weights_p), std::move(tails_p), std::move(heads_p));

        heads.clear();
        tails.clear();
        weights.clear();
        nodeWeights.clear();

        outputs[0] = G;
        outputs[1] = nodeWeightsArr;

    }
};

