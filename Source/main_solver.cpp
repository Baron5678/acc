#include "Graph.h"
#include "HungarianAlgorithm.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <string>
#include <iomanip>

using namespace std;

struct ApproxResult {
    int numCopies;          // Number of G copies found in H
    int totalExtEdges;      // Total edges added to extend H
    vector<vector<int>> extendedH;
    double hungarianTime;
};

struct SolveResult {
    Graph H_ext;
    double duration_sec = 0.0;

    int edgesG = 0;
    int edgesH = 0;
    int edgesHext = 0;

    int copiesRequested = 1;
    int copiesFound = 0;
    int totalEdgesAdded = 0;

    int bestDistance = INT_MAX;
    std::vector<int> bestMapping;
    bool isSubgraph = false;
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

SolveResult runHungarian(const Graph& G, const Graph& H, int targetCopies = -1) {
    SolveResult res{};
    res.edgesG = G.edgeCount();
    res.edgesH = H.edgeCount();

    res.copiesRequested = targetCopies;

    ApproxResult a = hungarianApproximateExtendMany(G, H, targetCopies);

    res.copiesFound = a.numCopies;
    res.totalEdgesAdded = a.totalExtEdges;

   
    res.duration_sec = a.hungarianTime / 1000.0;

    Graph H_ext(H.size);
    H_ext.adj = a.extendedH;

    res.H_ext = H_ext;
    res.edgesHext = H_ext.edgeCount();

    res.bestDistance = INT_MAX;
    res.isSubgraph = false;

    return res;
}

SolveResult runExact(const Graph& G, const Graph& H, int targetCopies) {
    using namespace std;
    using namespace std::chrono;

    SolveResult res{};
    res.edgesG = G.edgeCount();
    res.edgesH = H.edgeCount();

    auto start = high_resolution_clock::now();

    if (targetCopies == 1) {
        auto result = G.findBestMapping(H);
        res.bestMapping = result.first;
        res.bestDistance = result.second;

        Graph H_ext = H;
        if (res.bestDistance != 0 && res.bestDistance != INT_MAX) {
            for (int uG = 0; uG < G.size; ++uG) {
                for (int vG = 0; vG < G.size; ++vG) {
                    if (G.adj[uG][vG] > 0) {
                        int uH = res.bestMapping[uG];
                        int vH = res.bestMapping[vG];
                        if (H_ext.adj[uH][vH] == 0) {
                            H_ext.adj[uH][vH] = 1;
                        }
                    }
                }
            }
        }

        res.isSubgraph = (res.bestDistance == 0);
        res.copiesFound = (res.bestDistance == INT_MAX ? 0 : 1);
        res.H_ext = H_ext;
        res.edgesHext = H_ext.edgeCount();
        res.totalEdgesAdded = res.edgesHext - res.edgesH;

    }
    else {
        Graph H_ext = H;
        vector<bool> usedH(H.size, false);
        int totalEdgesAdded = 0;
        int copiesFound = 0;

        for (int copyNum = 0; copyNum < targetCopies; ++copyNum) {
            int availableNodes = 0;
            for (bool used : usedH) if (!used) availableNodes++;
            if (availableNodes < G.size) break;

            vector<int> availableH;
            for (int i = 0; i < H.size; ++i) {
                if (!usedH[i]) availableH.push_back(i);
            }

            Graph availableGraph((int)availableH.size());
            for (size_t i = 0; i < availableH.size(); ++i) {
                for (size_t j = 0; j < availableH.size(); ++j) {
                    availableGraph.adj[i][j] = H_ext.adj[availableH[i]][availableH[j]];
                }
            }

            auto result = G.findBestMapping(availableGraph);
            if (result.second == INT_MAX) break;

            vector<int> actualMapping(G.size);
            for (int i = 0; i < G.size; ++i) {
                actualMapping[i] = availableH[result.first[i]];
                usedH[actualMapping[i]] = true;
            }

            int edgesAdded = 0;
            for (int uG = 0; uG < G.size; ++uG) {
                for (int vG = 0; vG < G.size; ++vG) {
                    if (G.adj[uG][vG] > 0) {
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

        res.H_ext = H_ext;
        res.edgesHext = H_ext.edgeCount();
        res.totalEdgesAdded = totalEdgesAdded;
        res.copiesFound = copiesFound;

        res.isSubgraph = false;
        res.bestDistance = INT_MAX;
    }

    auto end = high_resolution_clock::now();
    res.duration_sec = duration_cast<duration<double>>(end - start).count();
    return res;
}

void displayResultsForSmallGraphs(std::string algo, Graph G, Graph H, Graph H_ext, double duration) {
    cout << "=== " << (algo == "exact" ? "EXACT" : "HUNGARIAN") << " ALGORITHM RESULTS ===" << endl;
    cout << "Algorithm time: " << fixed << setprecision(6) << duration << "ms" << endl << endl;

    cout << "-- Graph G --" << endl;
    G.print();
    cout << endl;
    cout << "-- H vs H extended --" << endl ;
    H_ext.printHighlighted(H);

    int addedEdges = 0;
    int n = std::max(H.size, H_ext.size);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int oldVal = (i < H.size && j < H.size) ? H.adj[i][j] : 0;
            int newVal = (i < H_ext.size && j < H_ext.size) ? H_ext.adj[i][j] : 0;
            if (newVal != 0 && oldVal == 0) addedEdges++;
        }
    }

    cout << "\nSummary:" << endl;
    cout << "Added edges (H_ext \\ H): " << addedEdges << endl;
    cout << "H size: " << H.edgeCount() << ", H_ext size: " << H_ext.edgeCount() << endl;

    cout << "==============================================" << endl;

}

void displayResultsForBigGraphs(std::string algo, int EdgesG, int EdgesH, int EdgesH_ext, double duration) {
    cout << "=== " << (algo == "exact" ? "EXACT" : "HUNGARIAN") << " ALGORITHM RESULTS (BIG GRAPHS) ===" << endl;
    cout << "Algorithm time: " << fixed << setprecision(6) << duration << " s" << endl;

    cout << "\nEdge counts:" << endl;
    cout << "  |E(G)|      = " << EdgesG << endl;
    cout << "  |E(H)|      = " << EdgesH << endl;
    cout << "  |E(H extended)|  = " << EdgesH_ext << endl;

    const int addedEdges = EdgesH_ext - EdgesH;
    cout << "\nExtension summary:" << endl;
    cout << "  Added edges (|E(H extended)| - |E(H)|) = " << addedEdges << endl;

    if (EdgesH > 0) {
        const double pct = 100.0 * (static_cast<double>(addedEdges) / static_cast<double>(EdgesH));
        cout << "  Relative increase vs H            = " << fixed << setprecision(2) << pct << "%" << endl;
    }
    else {
        cout << "  Relative increase vs H            = N/A (|E(H)| = 0)" << endl;
    }

    cout << "============================================================" << endl;
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

    SolveResult res;

    if (algorithm == "exact") {
        if (targetCopies <= 0) {
            res = runExact(G, H, 1); 
        } else {
            res = runExact(G, H, targetCopies);
        }
    } else if (algorithm == "hungarian") {
        res = runHungarian(G, H, targetCopies);
    } else {
        cerr << "Unknown algorithm: " << algorithm << endl;
        return 1;
    }

    if (G.size > 20) {
        displayResultsForBigGraphs(algorithm, res.edgesG, res.edgesH, res.edgesHext, res.duration_sec);

        if (algorithm == "hungarian" || (algorithm == "exact" && (targetCopies > 1))) {
            cout << "Copies requested: " << (res.copiesRequested <= 0 ? -1 : res.copiesRequested) << "\n";
            cout << "Copies found:     " << res.copiesFound << "\n";
            cout << "Total edges added: " << res.totalEdgesAdded << "\n";
        }
    }
    else {
        displayResultsForSmallGraphs(algorithm, G, H, res.H_ext, res.duration_sec);

        // Optional copies info also for small graphs
        if (algorithm == "hungarian" || (algorithm == "exact" && (targetCopies > 1))) {
            cout << "Copies requested: " << (res.copiesRequested <= 0 ? -1 : res.copiesRequested) << "\n";
            cout << "Copies found:     " << res.copiesFound << "\n";
            cout << "Total edges added: " << res.totalEdgesAdded << "\n";
        }
    }

    return 0;
}
