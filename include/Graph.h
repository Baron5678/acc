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
    void printWithHighlight(const std::vector<int>& mapping) const;
    int computeDistance(const Graph& other, const std::vector<int>& mapping) const;
    std::pair<std::vector<int>, int> findBestMapping(const Graph& target) const;
    Graph extendGraph(int target_size) const;
    void exactMinExtendGraph(const Graph& target);
    void printWithHighlightNewEdges(const Graph& originalH) const;

    
private:
    void loadFromFile(const std::string& filename, bool first_graph);
};