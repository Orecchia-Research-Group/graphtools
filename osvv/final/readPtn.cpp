/*
CPP MATLAB function: readPtn

USAGE: [(int64) partitions] = readPtn((char) ptnFilename);

PURPOSE: Read in a partition from a filename.
         On each line the corresponding node contains each community it belongs to.

NOTES:
   - ASSUMING C LONG TYPE IS 64 BITS

*/

#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"


std::map<int64_t, std::vector<int64_t> > readPtn(std::string ptnFilename) {
    std::ifstream ptnFile(ptnFilename);
    if (!ptnFile.is_open()) {
        std::cerr << "Failed to open ptnFilename " << ptnFilename <<". Check path and permissions.\n";
        return std::map<int64_t, std::vector<int64_t> >();
    }
    std::map<int64_t, std::vector<int64_t> > partitions;
    std::string line;
    int64_t node = 1;
    while(std::getline(ptnFile, line)) {
        std::istringstream iss{line};

        // Define variable 'data'. Use range constructor and stream iterator
        std::vector<int64_t> data{std::istream_iterator<int64_t>(iss), std::istream_iterator<int64_t>()};
        for (int64_t community: data) {
            partitions[community].push_back(node);
        }
        node++;
    }
    return partitions;
}

int main() {
    auto partitions = readPtn("/home/kameranis/Doctoral/graphtools/graphs/social_results/karate_final_100.ptn");
    for(auto part: partitions) {
        std::cout << "Partition " << part.first << ": " ;
        for(auto node: part.second) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }
}

class MexFunction : public matlab::mex::Function
{
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    matlab::data::ArrayFactory factory;

private:
    void mxErrorMessage(std::string msg)
    {
        matlab::data::CharArray args = factory.createCharArray(msg.data());
        matlabPtr->feval<void>(u"error", args);
    }

public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {
        std::string msg;
        if ((inputs.size() < 1) || (inputs.size() > 1)) {
            msg = "Error in readPtn. Wrong number of arguments.";
            mxErrorMessage(msg);
        }
        if ((outputs.size() < 1) || (outputs.size() > 1)) {
            msg = "Error in readPtn. Wrong number of return results.";
            mxErrorMessage(msg);
        }

        std::string ptnFilename(matlab::data::CharArray(inputs[0]).toAscii());
        auto temp = readPtn(ptnFilename);

        if (temp.size() == 0) {
            mxErrorMessage("Failed to read partitions.");
        }

        std::vector<matlab::data::Array> tempArray;
        for (auto part: temp) {
            size_t size = part.second.size();
            tempArray.push_back(factory.createArray({size, 1},
                    part.second.begin(), part.second.end()));
        }
        auto partitions = factory.createArray({1, temp.size()}, tempArray.begin(), tempArray.end());
        outputs[0] = partitions;
    }
};

