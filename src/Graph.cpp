#include "../include/Graph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <climits>
#include <functional>

using namespace std;

Graph::Graph() : size(0) {}

Graph::Graph(int n) : size(n) {
    adj.resize(n, vector<int>(n, 0));
}

Graph::Graph(const string& filename, bool first_graph) {
    loadFromFile(filename, first_graph);
}

void Graph::loadFromFile(const string& filename, bool first_graph) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        size = 0;
        return;
    }

    string line;
    int n_vertices = 0;

    // Skip first graph if reading second graph from file
    if (!first_graph) {
        if (getline(file, line)) {
            try {
                int n1 = stoi(line);
                // Skip adjacency matrix rows
                for (int i = 0; i < n1; ++i) {
                    getline(file, line);
                }
            }
            catch (...) {
                cerr << "Error parsing first graph size" << endl;
            }
        }
    }

    // Read target graph
    if (getline(file, line)) {
        try {
            n_vertices = stoi(line);
            size = n_vertices;
            adj.resize(size, vector<int>(size, 0));

            for (int i = 0; i < size; ++i) {
                if (getline(file, line)) {
                    stringstream ss(line);
                    int val;
                    for (int j = 0; j < size; ++j) {
                        if (ss >> val) {
                            adj[i][j] = val;
                        }
                    }
                }
            }
        }
        catch (...) {
            cerr << "Error parsing graph data" << endl;
            size = 0;
        }
    }
    file.close();
}

void Graph::resize(int new_size) {
    if (new_size > size) {
        for (auto& row : adj) {
            row.resize(new_size, 0);
        }
        adj.resize(new_size, vector<int>(new_size, 0));
        size = new_size;
    }
}

void Graph::printWithHighlight(const vector<int>& mapping) const {
    const string GREEN = "\033[1;32m";
    const string RESET = "\033[0m";

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout << adj[i][j] << " ";
        }
        cout << endl;
    }
}

void Graph::printWithHighlightNewEdges(const Graph& originalH) const {
    const string GREEN = "\033[1;32m";
    const string RESET = "\033[0m";

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            bool isNew = (adj[i][j] == 1 && originalH.adj[i][j] == 0);

            if (isNew) cout << GREEN;
            cout << adj[i][j];
            if (isNew) cout << RESET;

            cout << " ";
        }
        cout << "\n";
    }
}


int Graph::computeDistance(const Graph& other, const vector<int>& mapping) const {
    int cost = 0;
    int n = this->size;

    for (int uG = 0; uG < n; ++uG) {
        for (int vG = 0; vG < n; ++vG) {
            if (this->adj[uG][vG] > 0) {
                int uH = mapping[uG];
                int vH = mapping[vG];
                if (other.adj[uH][vH] == 0) {
                    cost++;
                }
            }
        }
    }
    return cost;
}

pair<vector<int>, int> Graph::findBestMapping(const Graph& target) const {
    vector<int> mapping(size, -1);
    vector<bool> usedH(target.size, false);
    
    // Recursive backtracking to find optimal mapping
    // Tries all possible assignments of G nodes to H nodes
    function<pair<vector<int>, int>(int)> solve = [&](int vG) -> pair<vector<int>, int> {
        if (vG == this->size) {
            // All G nodes mapped, compute extension cost
            int distance = this->computeDistance(target, mapping);
            return {mapping, distance};
        }

        int bestDistance = INT_MAX;
        vector<int> bestMapping;

        // Try mapping G node vG to each available H node
        for (int vH = 0; vH < target.size; ++vH) {
            if (!usedH[vH]) {
                mapping[vG] = vH;
                usedH[vH] = true;

                auto candidate = solve(vG + 1);
                if (candidate.second < bestDistance) {
                    bestDistance = candidate.second;
                    bestMapping = candidate.first;
                }

                // Backtrack
                usedH[vH] = false;
                mapping[vG] = -1;
            }
        }

        return {bestMapping, bestDistance};
    };
    
    return solve(0);
}

Graph Graph::extendGraph(int target_size) const {
    Graph extended(target_size);
    for (int i = 0; i < min(size, target_size); ++i) {
        for (int j = 0; j < min(size, target_size); ++j) {
            extended.adj[i][j] = adj[i][j];
        }
    }

    return extended;
}

void Graph::exactMinExtendGraph(const Graph& target) {
    Graph originalH = target;

    auto start = chrono::high_resolution_clock::now();

    auto result = findBestMapping(target);
    vector<int> bestMapping = result.first;
    int bestDistance = result.second;

    Graph H_ext = target;
    if (bestDistance != 0 && bestDistance != INT_MAX) {
        // Apply optimal mapping by adding missing edges to H
        for (int uG = 0; uG < this->size; ++uG) {
            for (int vG = 0; vG < this->size; ++vG) {
                if (this->adj[uG][vG] > 0) {
                    int uH = bestMapping[uG];
                    int vH = bestMapping[vG];
                    if (H_ext.adj[uH][vH] == 0) {
                        H_ext.adj[uH][vH] = 1;
                    }
                }
            }
        }
    }

    bool isSubgraph = (bestDistance == 0);

    auto end = chrono::high_resolution_clock::now();
    auto elapsed_ms = chrono::duration_cast<chrono::milliseconds>(end - start).count();

    cout << "Is Subgraph? " << (isSubgraph ? "YES" : "NO") << endl;
    cout << "Minimal Extension Cost: " << bestDistance << endl;

    if (bestDistance != INT_MAX) {
        cout << "One of the best mappings (G->H): ";
        for (size_t i = 0; i < bestMapping.size(); i++)
            cout << i << "->" << bestMapping[i] << " ";
        cout << endl;
    }

    cout << "Execution time: " << elapsed_ms << " ms" << endl;

    cout << endl << "--- Graph G ---" << endl;
    this->printWithHighlight(vector<int>());

    cout << endl << "--- Graph H ---" << endl;
    H_ext.printWithHighlightNewEdges(originalH);
}



