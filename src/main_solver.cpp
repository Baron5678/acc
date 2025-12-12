#include "../include/Graph.h"
#include "../include/HungarianAlgorithm.h"
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
    double greedyTime;
};

// Compute out-degree for each node in graph
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


pair<bool, vector<int>> hungarianMappingOne(const Graph& G, const Graph& H, const vector<bool>& usedH) {
    int n = G.size;
    int m = H.size;

    if (n > m) {
        return { false, {} };
    }

    // pack multiple disjoint copies, there must be at least n unused vertices.
    int unused = 0;
    for (bool u : usedH) if (!u) unused++;
    if (unused < n) {
        return { false, {} };
    }

    HungarianAlgorithm hungarian(m);

    // Build cost matrix based on structural similarity
    vector<int> degG = computeDegrees(G.adj);
    vector<int> degH = computeDegrees(H.adj);

    //Including dummy rows to make the cost matrix square, when m is different that n
    const int FORBIDDEN = 1'000'000;

    // Filling full m×m matrix: real rows use degree costs; dummy rows are all zeros.
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            if (i < n) {
                // Real row to forbid mapping into already-used vertices.
                if (usedH[j]) {
                    hungarian.setCost(i, j, FORBIDDEN);
                }
                else {
                    int cost = abs(degG[i] - degH[j]) + 1;
                    hungarian.setCost(i, j, cost);
                }
            }
            else {
                // Dummy row maps to 0 cost to any column.
                hungarian.setCost(i, j, 0);
            }
        }
    }

    vector<int> assignment = hungarian.findMinCostAssignment();
    vector<int> mapping(n);
    for (int i = 0; i < n; ++i) {
        mapping[i] = assignment[i];
        //ensure not hitting unused vertex.
        if (mapping[i] < 0 || mapping[i] >= m || usedH[mapping[i]]) {
            return { false, {} };
        }
    }

    return { true, mapping };
}

// Hungarian-based approximation: find multiple non-overlapping copies of G in H
ApproxResult hungarianApproximateExtendMany(const Graph& G, const Graph& H) {
    ApproxResult result;
    result.numCopies = 0;
    result.totalExtEdges = 0;
    result.extendedH = H.adj;
    
    auto start = chrono::high_resolution_clock::now();

    int m = H.size;
    vector<bool> usedH(m, false);
    
    while (true) {
        // Check if enough nodes remain for another copy
        int available = 0;
        for (bool used : usedH) {
            if (!used) available++;
        }
        
        if (available < G.size) break;

        // Find optimal mapping using current extended H
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

        // Mark nodes as used
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

// Greedy approximation: find multiple copies using simple degree matching
ApproxResult greedyApproximateExtendMany(const Graph& G, const Graph& H) {
    ApproxResult result;
    result.numCopies = 0;
    result.totalExtEdges = 0;
    result.extendedH = H.adj;
    
    auto start = chrono::high_resolution_clock::now();

    int m = H.size;
    vector<bool> usedH(m, false);
    vector<int> degG = computeDegrees(G.adj);

    while (true) {
        // Check if enough nodes remain for another copy
        int available = 0;
        for (bool used : usedH) {
            if (!used) available++;
        }
        
        if (available < G.size) break;

        vector<int> degH = computeDegrees(result.extendedH);
        vector<int> mapping(G.size, -1);
        vector<bool> tempUsed = usedH;
        
        // Greedy assignment: match nodes with most similar degrees
        for (int i = 0; i < G.size; ++i) {
            int bestH = -1;
            int bestScore = INT_MAX;

            for (int j = 0; j < m; ++j) {
                if (!tempUsed[j]) {
                    int score = abs(degG[i] - degH[j]);
                    if (score < bestScore) {
                        bestScore = score;
                        bestH = j;
                    }
                }
            }

            if (bestH == -1) break;

            mapping[i] = bestH;
            tempUsed[bestH] = true;
        }

        bool complete = true;
        for (int node : mapping) {
            if (node == -1) {
                complete = false;
                break;
            }
        }

        if (!complete) break;

        usedH = tempUsed;

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
    result.greedyTime = chrono::duration<double, milli>(end - start).count();

    return result;
}

void runApproximation(const Graph& G, const Graph& H) {
    cout << "\n=== HUNGARIAN ALGORITHM RESULTS ===" << endl;
    
    // Count original edge counts for comparison
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



    ApproxResult hungarianResult = hungarianApproximateExtendMany(G, H);

    cout << "Size of G (edges):          " << gEdges << endl;
    cout << "Size of H before (edges):   " << hEdges << endl;

    // Count edges in extended H
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
   
    cout << "G" << endl;

    G.printWithHighlightNewEdges(G);

    cout << endl;

    cout << "Original H" << endl;

    H.printWithHighlightNewEdges(H);

    cout << endl;

    Graph H_ext(H.size);
    H_ext.adj = hungarianResult.extendedH;

    cout << "\nExtended H adjacency matrix (new edges in green):\n";
    H_ext.printWithHighlightNewEdges(H);   // compare against original H


    // Run greedy for comparison
    ApproxResult greedyResult = greedyApproximateExtendMany(G, H);

    cout << "\n=== GREEDY COMPARISON ===" << endl;
    cout << "Greedy edges added:         " << greedyResult.totalExtEdges << endl;
    cout << "Greedy copies found:        " << greedyResult.numCopies << endl;
    cout << "Greedy time:                " << greedyResult.greedyTime << " ms" << endl;

    cout << "\n=== ALGORITHM COMPARISON ===" << endl;
    if (hungarianResult.totalExtEdges < greedyResult.totalExtEdges) {
        cout << "Hungarian vs Greedy edges:  " << hungarianResult.totalExtEdges
             << " vs " << greedyResult.totalExtEdges
             << " (Hungarian BETTER by " << (greedyResult.totalExtEdges - hungarianResult.totalExtEdges) << ")" << endl;
    } else if (hungarianResult.totalExtEdges > greedyResult.totalExtEdges) {
        cout << "Hungarian vs Greedy edges:  " << hungarianResult.totalExtEdges
             << " vs " << greedyResult.totalExtEdges
             << " (Greedy BETTER by " << (hungarianResult.totalExtEdges - greedyResult.totalExtEdges) << ")" << endl;
    } else {
        cout << "Hungarian vs Greedy edges:  " << hungarianResult.totalExtEdges
             << " vs " << greedyResult.totalExtEdges << " (TIE)" << endl;
    }

    cout << "Hungarian vs Greedy copies: " << hungarianResult.numCopies
         << " vs " << greedyResult.numCopies;
    if (hungarianResult.numCopies == greedyResult.numCopies) {
        cout << " (TIE)" << endl;
    } else {
        cout << endl;
    }

    cout << "Hungarian vs Greedy time:   " << hungarianResult.hungarianTime
         << " ms vs " << greedyResult.greedyTime << " ms" << endl;

}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " [algorithm] <input_file>" << endl;
        cerr << "Algorithms: exact | hungarian" << endl;
        return 1;
    }

    string algorithm = "hungarian";
    string inputFile;

    if (argc == 2) {
        inputFile = argv[1];
    } else {
        algorithm = argv[1];
        inputFile = argv[2];
    }

    transform(algorithm.begin(), algorithm.end(), algorithm.begin(), ::tolower);

    cout << "Algorithm: " << algorithm << endl;
    cout << "Input file: " << inputFile << endl;
    
    Graph G(inputFile, true);   // first graph
    Graph H(inputFile, false);  // second graph

    if (G.size == 0 || H.size == 0) {
        cerr << "Error: failed to load graphs from '" << inputFile << "'." << endl;
        return 1;
    }

    if (algorithm == "exact") {
        G.exactMinExtendGraph(H);
    } else if (algorithm == "hungarian") {
        runApproximation(G, H);
    } else {
        cerr << "Unknown algorithm: " << algorithm << endl;
        return 1;
    }

    return 0;
}
