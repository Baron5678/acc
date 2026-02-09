#pragma once

#include <vector>
#include <string>
#include <utility>

class Graph {
public:
    int size;
    std::vector<std::vector<int>> adj;

    Graph();
    Graph(int n);
    Graph(const std::string& filename, bool first_graph);
    
    void resize(int new_size);
    void print() const;
    void printHighlighted(const Graph& other) const;
    int ComputeDistance(const Graph& other, const std::vector<int>& mapping) const;
    int DistanceMetric(const Graph& other, const std::vector<int>& mapping) const;
    std::pair<std::vector<int>, int> FindBestMapping(const Graph& target) const;
    int edgeCount() const;
    static std::pair<bool, std::vector<int>> hungarianMappingOne(const Graph& G, const Graph& H);
    static std::pair<bool, std::vector<int>> hungarianMappingOne(const Graph& G, const Graph& H, const std::vector<bool>& usedH);
    
private:
    void loadFromFile(const std::string& filename, bool first_graph);
};