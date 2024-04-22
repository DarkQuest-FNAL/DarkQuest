#ifndef FUNCTIONS
#define FUNCTIONS

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "ROOT/RVec.hxx"

void BuildMap(std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> &trigMaps, std::string fname)
{
    /*
    read the trigger map from a txt file and store it in a 5D vector
    a[quadrant][dp11][dp12][dp21][dp22]
    4 quadrants, 80x80 for dp1, 50x50 for dp2
    */
    assert(trigMaps.size() == 4 && trigMaps[0].size() == 80 && trigMaps[0][0].size() == 80 && trigMaps[0][0][0].size() == 50 && trigMaps[0][0][0][0].size() == 50);

    // open the txt file and load the data
    std::ifstream file(fname);
    if (!file.is_open())
    {
        std::cout << "Error opening file: " << fname << std::endl;
        exit(1);
    }

    int nQuad = 0;
    std::string line;
    while (std::getline(file, line))
    {
        std::stringstream sline(line);
        std::vector<string> tokens;
        std::string token;
        while (sline >> token)
        {
            tokens.push_back(token);
        }
        if (tokens.size() > 0 && tokens[0] == "Quadrant")
        {
            std::cout << "Quadrant: " << tokens[1] << std::endl;
            nQuad = std::stoi(tokens[1]);
        }
        else if (tokens.size() > 0)
        {
            assert(tokens.size() == 5);
            int dp11 = std::stoi(tokens[0]);
            int dp12 = std::stoi(tokens[1]);
            int dp21 = std::stoi(tokens[2]);
            int dp22 = std::stoi(tokens[3]);
            float val = std::stof(tokens[4]);
            trigMaps[nQuad][dp11][dp12][dp21][dp22] = val;
        }
    }
}

bool IsTrigFired(std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> &trigMaps, const ROOT::VecOps::RVec<int> &hitIds, const ROOT::VecOps::RVec<int> &elmIds, const float threshold = 0.)
{
    /*
    check if the trigger is fired
    */
    std::vector<std::vector<int>> elmDPs(8, std::vector<int>());
    for (std::size_t i = 0; i < hitIds.size(); i++)
    {
        if (hitIds[i] < 55 || hitIds[i] > 62)
            continue;
        elmDPs[hitIds[i] - 55].push_back(elmIds[i] - 1);
    }
    for (std::size_t i = 0; i < 4; i++)
    {
        if (elmDPs[i].size() != 2 || elmDPs[i + 4].size() != 2)
            continue;
        std::sort(elmDPs[i].begin(), elmDPs[i].end());
        std::sort(elmDPs[i + 4].begin(), elmDPs[i + 4].end());
        int dp11 = elmDPs[i][0];
        int dp12 = elmDPs[i][1];
        int dp21 = elmDPs[i + 4][0];
        int dp22 = elmDPs[i + 4][1];
        if (trigMaps[i][dp11][dp12][dp21][dp22] > threshold)
            return 1;
    }
    return 0;
}

#endif
