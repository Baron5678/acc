#include "Graph.h"
#include "HungarianAlgorithm.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <string>

using namespace std;

struct ApproxResult {
    int numCopies;          // Number of G copies found in H
    int totalExtEdges;      // Total edges added to extend H
    vector<vector<int>> extendedH;
    double hungarianTime;
};

vector<int> computeDegrees(const vector<vector<int>>& adj);

int evaluateMapping(const Graph& G, const Graph& H, const vector<int>& mapping) { // Returns number of edges needed to add to H 
    //to accommodate G under the given mapping
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


pair<bool, vector<int>> hungarianMappingOne(const Graph& G, const Graph& H, const vector<bool>& usedH) {
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

    HungarianAlgorithm hungarian(m); //build square assignment matrix of size m x m

    vector<int> degG = computeDegrees(G.adj); //degrees of G nodes which are being mapped to H
    vector<int> degH = computeDegrees(H.adj); //degrees of H nodes
    const int FORBIDDEN = 1'000'000; //avoids overlapping copies
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


ApproxResult hungarianApproximateExtendMany(const Graph& G, const Graph& H, int targetCopies = -1) {
    ApproxResult result;
    result.numCopies = 0;
    result.totalExtEdges = 0;
    result.extendedH = H.adj;
    
    auto start = chrono::high_resolution_clock::now();

    int m = H.size;
    vector<bool> usedH(m, false);
    
    while (targetCopies == -1 || result.numCopies < targetCopies) {
        int available = 0;
        for (bool used : usedH) {
            if (!used) available++;
        }
        
        if (available < G.size) break;

        Graph tempH(m);
        tempH.adj = result.extendedH;
        
        auto mappingResult = hungarianMappingOne(G, tempH, usedH);
        if (!mappingResult.first) break;
        
        vector<int> mapping = mappingResult.second;

        bool validMapping = true;
        for (int node : mapping) {
            if (usedH[node]) {
                validMapping = false;
                break;
            }
        }

        if (!validMapping) break;

        for (int node : mapping) {
            usedH[node] = true;
        }

        int edgesAdded = 0;
        for (int i = 0; i < G.size; ++i) {
            for (int j = 0; j < G.size; ++j) {
                if (G.adj[i][j] == 1) {
                    int hi = mapping[i];
                    int hj = mapping[j];
                    if (result.extendedH[hi][hj] == 0) {
                        result.extendedH[hi][hj] = 1;
                        edgesAdded++;
                    }
                }
            }
        }

        result.totalExtEdges += edgesAdded;
        result.numCopies++;
    }

    auto end = chrono::high_resolution_clock::now();
    result.hungarianTime = chrono::duration<double, milli>(end - start).count();

    return result;
}

void runApproximation(const Graph& G, const Graph& H, int targetCopies = -1) {
    cout << "\n=== HUNGARIAN ALGORITHM RESULTS ===" << endl;
    if (targetCopies > 0) {
        cout << "Target copies requested: " << targetCopies << endl;
    } else {
        cout << "Target copies: maximum possible" << endl;
    }
    
    int gEdges = 0, hEdges = 0;
    for (int i = 0; i < G.size; ++i) {
        for (int j = 0; j < G.size; ++j) {
            gEdges += G.adj[i][j];
        }
    }
    for (int i = 0; i < H.size; ++i) {
        for (int j = 0; j < H.size; ++j) {
            hEdges += H.adj[i][j];
        }
    }

    ApproxResult hungarianResult = hungarianApproximateExtendMany(G, H, targetCopies);

    cout << "Size of G (edges):          " << gEdges << endl;
    cout << "Size of H before (edges):   " << hEdges << endl;


    int extendedEdges = 0;
    for (const auto& row : hungarianResult.extendedH) {
        for (int edge : row) {
            extendedEdges += edge;
        }
    }

    cout << "Size of H after (edges):    " << extendedEdges << endl;
    cout << "Edges added (extension size): " << hungarianResult.totalExtEdges << endl;
    cout << "Copies found:               " << hungarianResult.numCopies << endl;
    cout << "\nHungarian time: " << hungarianResult.hungarianTime << " ms" << endl;

    if (targetCopies > 0) {
        cout << "\n=== TARGET ACHIEVEMENT ===" << endl;
        cout << "Hungarian achieved target: " << (hungarianResult.numCopies >= targetCopies ? "YES" : "NO") 
             << " (" << hungarianResult.numCopies << "/" << targetCopies << ")" << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " [algorithm] <input_file> [number_of_copies]" << endl;
        cerr << "Algorithms: exact | hungarian" << endl;
        return 1;
    }

    string algorithm = "hungarian";
    string inputFile;
    int targetCopies = -1; // -1 means find maximum possible copies

    if (argc == 2) {
        inputFile = argv[1];
    } else if (argc == 3) {
        algorithm = argv[1];
        inputFile = argv[2];
    } else if (argc == 4) {
        algorithm = argv[1];
        inputFile = argv[2];
        targetCopies = stoi(argv[3]);
    }

    transform(algorithm.begin(), algorithm.end(), algorithm.begin(), ::tolower);

    cout << "Algorithm: " << algorithm << endl;
    cout << "Input file: " << inputFile << endl;
    if (targetCopies > 0) {
        cout << "Target copies: " << targetCopies << endl;
    } else {
        cout << "Target copies: maximum possible" << endl;
    }

    Graph G(inputFile, true);   // first graph
    Graph H(inputFile, false);  // second graph

    if (G.size == 0 || H.size == 0) {
        cerr << "Error: failed to load graphs from '" << inputFile << "'." << endl;
        return 1;
    }

    if (algorithm == "exact") {
        if (targetCopies <= 0) {
            G.exactMinExtendGraph(H, 1); // Default to single copy for exact algorithm
        } else {
            G.exactMinExtendGraph(H, targetCopies);
        }
    } else if (algorithm == "hungarian") {
        runApproximation(G, H, targetCopies);
    } else {
        cerr << "Unknown algorithm: " << algorithm << endl;
        return 1;
    }

    return 0;
}
