/* C MATLAB FUNCTION: loadmetisgraph

PURPOSE:    Reads in a METIS file as described in the manual 
            (http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf)
            and return the graph and the node weights as specified in the file

USAGE:      function [G, weights] = loadmetsgraph(graphFilename);

INPUTS:
-graphFilename (char):  Path to the graph file

Outputs:
-G (sparse matrix):     A sparse representation of the graph read

TODO: Be able to read all METIS formats
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


void loadMetisGraph(std::string graphFilename, size_t &n, size_t &m, std::vector<size_t> &heads, std::vector<size_t> &tails, std::vector<double> &weights, std::vector<double> &nodeWeights) {
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
    while (std::getline(graphFile, line)) {
        std::istringstream iss{line};
        // Define variable 'data'. Use range constructor and stream iterator
        std::vector<int64_t> data{std::istream_iterator<size_t>(iss), std::istream_iterator<size_t>()};
        
        if (n == 0) {
            n = data[0];
            m = data[1];
            if (data.size() > 2)
                flag = data[2];
            continue;
        }
        nodeWeights.push_back(data[0]);
        for (size_t j = 1; j < data.size(); j += 2) {
            heads.push_back(data[j] - 1);
            tails.push_back(node);
            weights.push_back((double) (data[j + 1]));
        }
        node++;
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
        loadMetisGraph(graphFilename, n, m, heads, tails, weights, nodeWeights);

        auto nodeWeightsArr = factory.createArray({1, n}, nodeWeights.begin(), nodeWeights.end());

        matlab::data::buffer_ptr_t<size_t> heads_p = factory.createBuffer<size_t>(2 * m);
        matlab::data::buffer_ptr_t<size_t> tails_p = factory.createBuffer<size_t>(2 * m);
        matlab::data::buffer_ptr_t<double> weights_p = factory.createBuffer<double>(2 * m);

        size_t *headsPtr = heads_p.get();
        size_t *tailsPtr = tails_p.get();
        double *weightsPtr = weights_p.get();


        std::for_each(heads.begin(), heads.end(), [&](const size_t& e) { *(headsPtr++) = e; });
        std::for_each(tails.begin(), tails.end(), [&](const size_t& e) { *(tailsPtr++) = e; });
        std::for_each(weights.begin(), weights.end(), [&](const double& e) { *(weightsPtr++) = e; });

        // Use the buffers to create the sparse array
        matlab::data::SparseArray<double> G =
            factory.createSparseArray<double>({n, n}, 2 * m,
                std::move(weights_p), std::move(tails_p), std::move(heads_p));

        heads.clear();
        tails.clear();
        weights.clear();
        nodeWeights.clear();

        outputs[0] = G;
        outputs[1] = nodeWeightsArr;

    }
};
