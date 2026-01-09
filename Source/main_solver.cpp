#include "Graph.h"
#include "HungarianAlgorithm.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <string>
#include <iomanip>
#include <climits>
#include <unordered_map>
#include <sstream>
#include <functional>

using namespace std;

struct BestPerSet {
    int dist;
    vector<int> mapping;
};

struct ApproxResult {
    int numCopies;         
    int totalExtEdges;     
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

struct Candidate {
    int dist;
    vector<int> mapping;
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


pair<bool, vector<int>> hungarianMappingOne(
    const Graph& G,
    const Graph& H,
    const vector<bool>* forbidColsRow0 = nullptr // if provided, columns marked true are forbidden for row 0
) {
    int n = G.size;
    int m = H.size;

    if (n > m) {
        return { false, {} };
    }

    HungarianAlgorithm hungarian(m); // build square assignment matrix of size m x m

    vector<int> degG = computeDegrees(G.adj);
    vector<int> degH = computeDegrees(H.adj);

    const int FORBIDDEN = 1'000'000;

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            if (i < n) {
                // If we want to force a different vertex-set than some previous one:
                // forbid ALL vertices of that set for row 0, guaranteeing ≥1 vertex differs.
                if (i == 0 && forbidColsRow0 && j < (int)forbidColsRow0->size() && (*forbidColsRow0)[j]) {
                    hungarian.setCost(i, j, FORBIDDEN);
                    continue;
                }

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
            else {
                hungarian.setCost(i, j, 0);
            }
        }
    }

    vector<int> assignment = hungarian.findMinCostAssignment();

    vector<int> mapping(n);
    for (int i = 0; i < n; ++i) {
        mapping[i] = assignment[i];
        if (mapping[i] < 0 || mapping[i] >= m) {
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

    const int n = G.size;
    const int m = H.size;

    if (n > m) {
        auto end = chrono::high_resolution_clock::now();
        result.hungarianTime = chrono::duration<double, milli>(end - start).count();
        return result;
    }

    // Store only the *vertex sets* used by accepted copies (order-independent).
    vector<vector<int>> previousVertexSets;

    while (targetCopies == -1 || result.numCopies < targetCopies) {
        Graph tempH(m);
        tempH.adj = result.extendedH;

        // 1) Get a Hungarian mapping
        auto mappingResult = hungarianMappingOne(G, tempH, nullptr);
        if (!mappingResult.first) break;

        vector<int> mapping = mappingResult.second;

        // 2) Normalize to a vertex-set (sorted) for the "distinct copy" rule
        vector<int> vertexSet = mapping;
        sort(vertexSet.begin(), vertexSet.end());

        // 3) If this vertex-set equals any previous one, force a change:
        //    forbid ALL vertices of that set for row 0 -> guarantees at least one vertex differs.
        int safety = 0;
        while (true) {
            bool duplicate = false;
            for (const auto& prevSet : previousVertexSets) {
                if (prevSet == vertexSet) { duplicate = true; break; }
            }
            if (!duplicate) break;

            vector<bool> forbidCols(m, false);
            for (int v : vertexSet) forbidCols[v] = true;

            auto altRes = hungarianMappingOne(G, tempH, &forbidCols);
            if (!altRes.first) {
                // No alternative mapping that differs by ≥1 vertex exists
                safety = 1000;
                break;
            }

            mapping = altRes.second;
            vertexSet = mapping;
            sort(vertexSet.begin(), vertexSet.end());

            if (++safety > 1000) break;
        }

        if (safety >= 1000) break;

        // 4) Accept mapping and extend H (edges may overlap; count only newly added edges)
        int edgesAdded = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
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

        previousVertexSets.push_back(vertexSet);
    }

    auto end = chrono::high_resolution_clock::now();
    result.hungarianTime = chrono::duration<double, milli>(end - start).count();

    return result;
}

void runApproximation(const Graph& G, const Graph& H, int targetCopies = -1) {
    cout << "\n=== HUNGARIAN ALGORITHM RESULTS ===" << endl;
    if (targetCopies > 0) {
        cout << "Target copies requested: " << targetCopies << endl;
    }
    else {
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

static int ExtendGraph(const Graph& G, Graph& H_ext, const std::vector<int>& mapping) {
    const int n = G.size;
    int added = 0;

    for (int uG = 0; uG < n; ++uG) {
        int uH = mapping[uG];
        for (int vG = 0; vG < n; ++vG) {
            if (G.adj[uG][vG] > 0) {
                int vH = mapping[vG];
                if (H_ext.adj[uH][vH] == 0) {
                    H_ext.adj[uH][vH] = 1;
                    added++;
                }
            }
        }
    }
    return added;
}


SolveResult ExactMinExtendGraph(const Graph& G, const Graph& H, int targetCopies) {
    using namespace std;
    using namespace std::chrono;

    SolveResult res{};
    res.edgesG = G.edgeCount();
    res.edgesH = H.edgeCount();

    auto start = high_resolution_clock::now();
    Graph H_ext = H;
    if (targetCopies == 1) {
        auto result = G.FindBestMapping(H);
        res.bestMapping = result.first;
        res.bestDistance = result.second;

        
        if (res.bestDistance != 0 && res.bestDistance != INT_MAX) {
            ExtendGraph(G, H_ext, res.bestMapping);
        }

        res.isSubgraph = (res.bestDistance == 0);
        res.copiesFound = (res.bestDistance == INT_MAX ? 0 : 1);
        res.H_ext = H_ext;
        res.edgesHext = H_ext.edgeCount();
        res.totalEdgesAdded = res.edgesHext - res.edgesH;

    }
    else {

        const int n = G.size;
        const int m = H.size;
        auto distanceUnderMapping = [&](const Graph& Hcur, const vector<int>& mapping) -> int {
            int dist = 0;
            for (int uG = 0; uG < n; ++uG) {
                int uH = mapping[uG];
                for (int vG = 0; vG < n; ++vG) {
                    if (G.adj[uG][vG] > 0) {
                        int vH = mapping[vG];
                        if (Hcur.adj[uH][vH] == 0) dist++;
                    }
                }
            }
            return dist;
            };

        auto vertexSetKey = [&](const vector<int>& mapping) -> string {
            vector<int> s = mapping;
            sort(s.begin(), s.end());
            ostringstream oss;
            for (int i = 0; i < (int)s.size(); ++i) {
                if (i) oss << ',';
                oss << s[i];
            }
            return oss.str();
            };

        

        unordered_map<string, BestPerSet> bestForSet;
        bestForSet.reserve(1024);
        vector<int> mapping(n, -1);
        vector<char> usedH(m, false);

        function<void(int)> dfs = [&](int uG) {
            if (uG == n) {
                int d = distanceUnderMapping(H, mapping);
                string key = vertexSetKey(mapping);

                auto it = bestForSet.find(key);
                if (it == bestForSet.end() || d < it->second.dist) {
                    bestForSet[key] = BestPerSet{ d, mapping };
                }
                return;
            }

            for (int vH = 0; vH < m; ++vH) {
                if (usedH[vH]) continue;
                usedH[vH] = true;
                mapping[uG] = vH;
                dfs(uG + 1);
                mapping[uG] = -1;
                usedH[vH] = false;
            }
            };

        dfs(0);

        vector<Candidate> cand;
        cand.reserve(bestForSet.size());
        for (auto& kv : bestForSet) {
            cand.push_back(Candidate{ kv.second.dist, kv.second.mapping });
        }

        sort(cand.begin(), cand.end(), [](const Candidate& a, const Candidate& b) {
            return a.dist < b.dist;
            });

        int copiesFound = min(targetCopies, (int)cand.size());

        int edgesAddedTotal = 0;

        for (int idx = 0; idx < copiesFound; ++idx) {
            const auto& mapGtoH = cand[idx].mapping;

            for (int uG = 0; uG < n; ++uG) {
                int uH = mapGtoH[uG];
                for (int vG = 0; vG < n; ++vG) {
                    if (G.adj[uG][vG] > 0) {
                        int vH = mapGtoH[vG];
                        if (H_ext.adj[uH][vH] == 0) {
                            H_ext.adj[uH][vH] = 1;
                            edgesAddedTotal++;
                        }
                    }
                }
            }
        }

        res.H_ext = H_ext;
        res.edgesHext = H_ext.edgeCount();
        res.totalEdgesAdded = edgesAddedTotal;
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
    cout << "-- H vs H extended --" << endl;
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
    }
    else if (argc == 3) {
        algorithm = argv[1];
        inputFile = argv[2];
    }
    else if (argc == 4) {
        algorithm = argv[1];
        inputFile = argv[2];
        targetCopies = stoi(argv[3]);
    }

    transform(algorithm.begin(), algorithm.end(), algorithm.begin(), ::tolower);

    cout << "Algorithm: " << algorithm << endl;
    cout << "Input file: " << inputFile << endl;
    if (targetCopies > 0) {
        cout << "Target copies: " << targetCopies << endl;
    }
    else {
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
            res = ExactMinExtendGraph(G, H, 1);
        }
        else {
            res = ExactMinExtendGraph(G, H, targetCopies);
        }
    }
    else if (algorithm == "hungarian") {
        res = runHungarian(G, H, targetCopies);
    }
    else {
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
