#pragma once

#include <vector>
#include <string>

class GraphGenerator {
public:
    static std::vector<std::vector<int>> generateConnectedGraph(int n, double density = 0.8);
    static void saveGraphsToFile(const std::string& filename,
                                const std::vector<std::vector<int>>& G,
                                const std::vector<std::vector<int>>& H);
};