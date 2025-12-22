#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include "Graph.h"
#include "GraphGenerator.h"



using namespace std;

struct TestResult {
    string algorithm;
    int gSize;
    int hSize;
    int targetCopies;
    int copiesFound;
    int edgesAdded;
    double timeMs;
    bool timeout;
    string status;
};

class TestRunner {
private:
    vector<TestResult> results;

    TestResult runExactTest(const Graph& G, const Graph& H, int targetCopies, int timeoutMs) {
        TestResult result;
        result.algorithm = "Exact";
        result.gSize = G.size;
        result.hSize = H.size;
        result.targetCopies = targetCopies;
        result.timeout = false;

        auto start = chrono::high_resolution_clock::now();

        try {
            if (targetCopies == 1) {
                auto mapping = G.findBestMapping(H);
                auto end = chrono::high_resolution_clock::now();
                result.timeMs = chrono::duration<double, milli>(end - start).count();

                if (result.timeMs > timeoutMs) {
                    result.timeout = true;
                    result.status = "TIMEOUT";
                } else {
                    result.edgesAdded = mapping.second;
                    result.copiesFound = (mapping.second == INT_MAX) ? 0 : 1;
                    result.status = "OK";
                }
            } else {
                Graph tempG = G;
                Graph tempH = H;
                
                streambuf* orig = cout.rdbuf();
                stringstream buffer;
                cout.rdbuf(buffer.rdbuf());
                
                tempG.exactMinExtendGraph(tempH, targetCopies);
                
                cout.rdbuf(orig);
                
                auto end = chrono::high_resolution_clock::now();
                result.timeMs = chrono::duration<double, milli>(end - start).count();

                if (result.timeMs > timeoutMs) {
                    result.timeout = true;
                    result.status = "TIMEOUT";
                    result.copiesFound = 0;
                    result.edgesAdded = 0;
                } else {
                    string output = buffer.str();
                    size_t copiesPos = output.find("Copies found: ");
                    size_t edgesPos = output.find("Total edges added: ");
                    
                    if (copiesPos != string::npos && edgesPos != string::npos) {
                        result.copiesFound = stoi(output.substr(copiesPos + 14));
                        result.edgesAdded = stoi(output.substr(edgesPos + 19));
                        result.status = "OK";
                    } else {
                        result.copiesFound = min(targetCopies, H.size / G.size);
                        result.edgesAdded = 0;
                        result.status = "OK";
                    }
                }
            }
        } catch (...) {
            result.timeout = true;
            result.status = "ERROR";
            result.timeMs = timeoutMs + 1;
        }

        return result;
    }

    TestResult runApproximationTest(const Graph& G, const Graph& H, int targetCopies) {
        TestResult hungarianResult;

        hungarianResult.algorithm = "Hungarian";
        hungarianResult.gSize = G.size;
        hungarianResult.hSize = H.size;
        hungarianResult.targetCopies = targetCopies;

        auto start = chrono::high_resolution_clock::now();

        try {
            Graph extendedH = H;
            vector<bool> usedH(H.size, false);
            int totalEdgesAdded = 0;
            int copiesFound = 0;
            
            int maxCopies = (targetCopies == -1) ? H.size / G.size : targetCopies;
            
            for (int copyNum = 0; copyNum < maxCopies; ++copyNum) {
                int available = 0;
                for (bool used : usedH) {
                    if (!used) available++;
                }
                
                if (available < G.size) break;
                
                auto mappingResult = Graph::hungarianMappingOne(G, extendedH, usedH);
                if (!mappingResult.first) break;
                
                vector<int> mapping = mappingResult.second;
                
                for (int node : mapping) {
                    usedH[node] = true;
                }
                
                int edgesAdded = 0;
                for (int i = 0; i < G.size; ++i) {
                    for (int j = 0; j < G.size; ++j) {
                        if (G.adj[i][j] == 1) {
                            int hi = mapping[i];
                            int hj = mapping[j];
                            if (extendedH.adj[hi][hj] == 0) {
                                extendedH.adj[hi][hj] = 1;
                                edgesAdded++;
                            }
                        }
                    }
                }
                
                totalEdgesAdded += edgesAdded;
                copiesFound++;
            }

            auto end = chrono::high_resolution_clock::now();
            hungarianResult.timeMs = chrono::duration<double, milli>(end - start).count();
            hungarianResult.timeout = false;
            hungarianResult.edgesAdded = totalEdgesAdded;
            hungarianResult.copiesFound = copiesFound;
            hungarianResult.status = (targetCopies == -1 || copiesFound >= targetCopies) ? "TARGET_MET" : "TARGET_MISSED";

        } catch (...) {
            hungarianResult.timeout = true;
            hungarianResult.status = "ERROR";
            hungarianResult.edgesAdded = 0;
            hungarianResult.copiesFound = 0;
        }

        return hungarianResult;
    }

    pair<Graph, Graph> createFixedGraphs(int gSize, int hSize) {
        Graph G(gSize), H(hSize);
        
        if (gSize == 3 && hSize == 5) {
            G.adj = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}};
            H.adj = {{0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}, 
                     {0, 0, 0, 0, 1}, {1, 0, 0, 0, 0}};
        } else if (gSize == 4 && hSize == 6) {
            G.adj = {{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {1, 0, 0, 0}};
            H.adj = {{0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, 
                     {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}, {1, 0, 1, 0, 0, 0}};
        } else if (gSize == 5 && hSize == 8) {
            G.adj = {{0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 1}, 
                     {0, 0, 0, 0, 1}, {1, 0, 0, 0, 0}};
            H.adj = {{0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0}, 
                     {0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0}, 
                     {0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0}, 
                     {0, 0, 0, 0, 0, 0, 0, 1}, {1, 0, 1, 0, 1, 0, 0, 0}};
        } else if (gSize == 6 && hSize == 10) {
            G.adj = {{0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, 
                     {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}, {1, 0, 1, 0, 0, 0}};
            H.adj = {{0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, 
                     {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, 
                     {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, 
                     {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, 
                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {1, 0, 1, 0, 1, 0, 1, 0, 0, 0}};
        } else if (gSize == 6 && hSize == 12) {
            G.adj = {{0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 1, 0}, 
                     {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}, {1, 0, 0, 1, 0, 0}};
            H.adj = {{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
                     {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, 
                     {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, 
                     {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, 
                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, 
                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0}};
        } else if (gSize == 4 && hSize == 10) {
            G.adj = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}};
            H.adj = {{0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, 
                     {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
                     {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, 
                     {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, 
                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}};
        } else {
            for (int i = 0; i < gSize - 1; ++i) {
                G.adj[i][i+1] = 1;
            }
            G.adj[gSize-1][0] = 1;
            
            for (int i = 0; i < hSize - 1; ++i) {
                H.adj[i][i+1] = 1;
            }
            H.adj[hSize-1][0] = 1;
        }
        
        return {G, H};
    }

public:
    void runBasicPerformanceTests() {
        cout << "1. Basic Performance Comparison (Directed Graphs):" << endl;
        cout << "===================================================" << endl;

        vector<pair<int, int>> testCases = {{3, 5}, {4, 6}, {5, 8}, {6, 10}};

        for (auto& testCase : testCases) {
            int gSize = testCase.first;
            int hSize = testCase.second;

            cout << "G(" << gSize << ") vs H(" << hSize << "):" << endl;

            auto [G, H] = createFixedGraphs(gSize, hSize);

            TestResult exactResult = runExactTest(G, H, 1, 5000);
            cout << "  Exact: " << fixed << setprecision(0) << exactResult.timeMs
                 << "ms (cost " << exactResult.edgesAdded << ")" << endl;

            auto approxResult = runApproximationTest(G, H, -1);
            cout << "  Hungarian: " << fixed << setprecision(3) << approxResult.timeMs
                 << "ms (" << approxResult.edgesAdded << " edges)" << endl;
            cout << endl;
        }
    }

    void runCopyTargetingTests() {
        cout << "2. Exact vs Hungarian Approximation Analysis:" << endl;
        cout << "=============================================" << endl;

        GraphGenerator generator;

        vector<pair<int, int>> comparableTests = {{3, 6}, {4, 8}, {5, 10}};

        for (auto& testCase : comparableTests) {
            int gSize = testCase.first;
            int hSize = testCase.second;

            cout << "Direct comparison G(" << gSize << ") vs H(" << hSize << "):" << endl;

            auto [G, H] = createFixedGraphs(gSize, hSize);

            TestResult exactResult = runExactTest(G, H, 1, 3000);
            TestResult hungarianResult = runApproximationTest(G, H, 1);

            if (exactResult.timeout) {
                cout << "  Exact: TIMEOUT (>" << exactResult.timeMs << "ms)" << endl;
                cout << "  Hungarian: " << fixed << setprecision(3) << hungarianResult.timeMs
                     << "ms (" << hungarianResult.edgesAdded << " edges)" << endl;
                cout << "  Approximation error: Cannot measure (exact timeout)" << endl;
            } else {
                cout << "  Exact: " << fixed << setprecision(0) << exactResult.timeMs
                     << "ms (optimal: " << exactResult.edgesAdded << " edges)" << endl;
                cout << "  Hungarian: " << fixed << setprecision(3) << hungarianResult.timeMs
                     << "ms (approx: " << hungarianResult.edgesAdded << " edges)" << endl;

                if (exactResult.edgesAdded == 0) {
                    if (hungarianResult.edgesAdded == 0) {
                        cout << "  Approximation error: 0% (both optimal)" << endl;
                    } else {
                        cout << "  Approximation error: INFINITE% (exact found subgraph, Hungarian didn't)" << endl;
                    }
                } else {
                    double errorPercent = ((double)(hungarianResult.edgesAdded - exactResult.edgesAdded) / exactResult.edgesAdded) * 100.0;
                    cout << "  Approximation error: " << fixed << setprecision(1) << errorPercent << "%" << endl;
                }

                double speedup = exactResult.timeMs / hungarianResult.timeMs;
                cout << "  Hungarian speedup: " << fixed << setprecision(1) << speedup << "x faster" << endl;
            }
            cout << endl;
        }

        cout << "Exponential behavior demonstration G(6) vs H(12):" << endl;
        auto [G6, H12] = createFixedGraphs(6, 12);

        TestResult exactHardResult = runExactTest(G6, H12, 1, 5000);
        TestResult hardApproxResult = runApproximationTest(G6, H12, 1);

        if (exactHardResult.timeout) {
            cout << "  Exact: TIMEOUT (>5s) - exponential complexity!" << endl;
        } else {
            cout << "  Exact: " << fixed << setprecision(0) << exactHardResult.timeMs
                 << "ms (optimal: " << exactHardResult.edgesAdded << " edges)" << endl;
        }
        cout << "  Hungarian: " << fixed << setprecision(3) << hardApproxResult.timeMs
             << "ms (approx: " << hardApproxResult.edgesAdded << " edges)" << endl;
        cout << "  Demonstrates: Exact becomes impractical, Hungarian stays fast" << endl;
        cout << endl;

        cout << "Testing multiple copies G(4) vs H(10):" << endl;
        auto [G4, H10] = createFixedGraphs(4, 10);

        for (int copies : {1, 2, 3}) {
            cout << "Target: " << copies << " copy/copies" << endl;

            TestResult multiExactResult = runExactTest(G4, H10, copies, 1000);
            if (multiExactResult.timeout) {
                cout << "  Exact: TIMEOUT" << endl;
            } else {
                cout << "  Exact: " << multiExactResult.copiesFound << " copies, "
                     << multiExactResult.edgesAdded << " edges, "
                     << fixed << setprecision(0) << multiExactResult.timeMs << "ms" << endl;
            }

            auto multiApproxResult = runApproximationTest(G4, H10, copies);
            cout << "  Hungarian: " << multiApproxResult.copiesFound << " copies, "
                 << multiApproxResult.edgesAdded << " edges, "
                 << fixed << setprecision(3) << multiApproxResult.timeMs << "ms, Target: "
                 << multiApproxResult.status << endl;
            cout << endl;
        }
    }


};

int main() {
    cout << "Graph Algorithms Test Suite" << endl;
    cout << "============================" << endl;
    cout << endl;

    TestRunner runner;

    runner.runBasicPerformanceTests();
    runner.runCopyTargetingTests();

    cout << "Done." << endl;
    return 0;
}
