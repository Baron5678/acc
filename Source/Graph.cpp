#include "Graph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <climits>
#include <functional>
#include "HungarianAlgorithm.h"

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
    
    function<pair<vector<int>, int>(int)> solve = [&](int vG) -> pair<vector<int>, int> {
        if (vG == this->size) {
            int distance = this->computeDistance(target, mapping);
            return {mapping, distance};
        }

        int bestDistance = INT_MAX;
        vector<int> bestMapping;

        for (int vH = 0; vH < target.size; ++vH) {
            if (!usedH[vH]) {
                mapping[vG] = vH;
                usedH[vH] = true;

                auto candidate = solve(vG + 1);
                if (candidate.second < bestDistance) {
                    bestDistance = candidate.second;
                    bestMapping = candidate.first;
                }

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

void Graph::exactMinExtendGraph(const Graph& target, int targetCopies) {
    Graph originalH = target;

    auto start = chrono::high_resolution_clock::now();

    if (targetCopies == 1) {
        auto result = findBestMapping(target);
        vector<int> bestMapping = result.first;
        int bestDistance = result.second;

        Graph H_ext = target;
        if (bestDistance != 0 && bestDistance != INT_MAX) {
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
        H_ext.printWithHighlight(vector<int>());
    } else {
        Graph H_ext = target;
        vector<bool> usedH(target.size, false);
        int totalEdgesAdded = 0;
        int copiesFound = 0;
        
        for (int copyNum = 0; copyNum < targetCopies; ++copyNum) {
            int availableNodes = 0;
            for (bool used : usedH) {
                if (!used) availableNodes++;
            }
            
            if (availableNodes < this->size) {
                break;
            }
            
            vector<int> availableH;
            for (int i = 0; i < target.size; ++i) {
                if (!usedH[i]) {
                    availableH.push_back(i);
                }
            }
            
            Graph availableGraph(availableH.size());
            for (size_t i = 0; i < availableH.size(); ++i) {
                for (size_t j = 0; j < availableH.size(); ++j) {
                    availableGraph.adj[i][j] = H_ext.adj[availableH[i]][availableH[j]];
                }
            }
            
            auto result = findBestMapping(availableGraph);
            if (result.second == INT_MAX) {
                break;
            }
            
            vector<int> actualMapping(this->size);
            for (int i = 0; i < this->size; ++i) {
                actualMapping[i] = availableH[result.first[i]];
                usedH[actualMapping[i]] = true;
            }
            
            int edgesAdded = 0;
            for (int uG = 0; uG < this->size; ++uG) {
                for (int vG = 0; vG < this->size; ++vG) {
                    if (this->adj[uG][vG] > 0) {
                        int uH = actualMapping[uG];
                        int vH = actualMapping[vG];
                        if (H_ext.adj[uH][vH] == 0) {
                            H_ext.adj[uH][vH] = 1;
                            edgesAdded++;
                        }
                    }
                }
            }
            
            totalEdgesAdded += edgesAdded;
            copiesFound++;
        }

        auto end = chrono::high_resolution_clock::now();
        auto elapsed_ms = chrono::duration_cast<chrono::milliseconds>(end - start).count();

        cout << "\n=== EXACT ALGORITHM RESULTS (MULTIPLE COPIES) ===" << endl;
        cout << "Target copies requested: " << targetCopies << endl;
        cout << "Copies found: " << copiesFound << endl;
        cout << "Total edges added: " << totalEdgesAdded << endl;
        cout << "Execution time: " << elapsed_ms << " ms" << endl;

        cout << endl << "--- Graph G ---" << endl;
        this->printWithHighlight(vector<int>());

        cout << endl << "--- Extended Graph H ---" << endl;
        H_ext.printWithHighlight(vector<int>());
    }
}

vector<int> computeDegrees(const vector<vector<int>>& adj) {
    int n = adj.size();
    vector<int> degrees(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            degrees[i] += adj[i][j];
        }
    }
    return degrees;
}

static int evaluateMapping(const Graph& G, const Graph& H, const vector<int>& mapping) {
    int edgesNeeded = 0;
    for (int i = 0; i < G.size; ++i) {
        for (int j = 0; j < G.size; ++j) {
            if (G.adj[i][j] == 1) {
                int hi = mapping[i];
                int hj = mapping[j];
                if (H.adj[hi][hj] == 0) {
                    edgesNeeded++;
                }
            }
        }
    }
    return edgesNeeded;
}

pair<bool, vector<int>> Graph::hungarianMappingOne(const Graph& G, const Graph& H) {
    vector<bool> usedH(H.size, false);
    return hungarianMappingOne(G, H, usedH);
}

pair<bool, vector<int>> Graph::hungarianMappingOne(const Graph& G, const Graph& H, const vector<bool>& usedH) {
    int n = G.size;
    int m = H.size;

    if (n > m) {
        return { false, {} };
    }

    int unused = 0;
    for (bool u : usedH) if (!u) unused++;
    if (unused < n) {
        return { false, {} };
    }

    HungarianAlgorithm hungarian(m);

    vector<int> degG = computeDegrees(G.adj);
    vector<int> degH = computeDegrees(H.adj);
    const int FORBIDDEN = 1'000'000;

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            if (i < n) {
                if (usedH[j]) {
                    hungarian.setCost(i, j, FORBIDDEN);
                }
                else {
                    int cost = 0;
                    
                    for (int k = 0; k < n; ++k) {
                        if (i != k) {
                            if (G.adj[i][k] == 1) {
                                int outDegreeH = 0;
                                for (int l = 0; l < m; ++l) {
                                    if (H.adj[j][l] == 1) outDegreeH++;
                                }
                                if (outDegreeH == 0) cost += 20;
                            }
                            if (G.adj[k][i] == 1) {
                                int inDegreeH = 0;
                                for (int l = 0; l < m; ++l) {
                                    if (H.adj[l][j] == 1) inDegreeH++;
                                }
                                if (inDegreeH == 0) cost += 20;
                            }
                        }
                    }
                    
                    cost += abs(degG[i] - degH[j]);
                    cost += (i + j) / 10;
                    
                    hungarian.setCost(i, j, cost + 1);
                }
            }
            else {
                hungarian.setCost(i, j, 0);
            }
        }
    }

    vector<int> assignment = hungarian.findMinCostAssignment();
    vector<int> mapping(n);
    for (int i = 0; i < n; ++i) {
        mapping[i] = assignment[i];
        if (mapping[i] < 0 || mapping[i] >= m || usedH[mapping[i]]) {
            return { false, {} };
        }
    }

    return { true, mapping };
}
