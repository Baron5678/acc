#pragma once

#include <vector>

class HungarianAlgorithm {
public:
    HungarianAlgorithm(int size);
    
    void setCost(int i, int j, int cost);
    std::vector<int> findMinCostAssignment();
    
private:
    std::vector<std::vector<int>> cost_matrix;
    int n;
    
    std::vector<int> solve();
};