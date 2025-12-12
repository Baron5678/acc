#include "../include/GraphGenerator.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <numeric>

using namespace std;

vector<vector<int>> GraphGenerator::generateConnectedGraph(int n, double density) {
    if (n <= 0) return {};

    vector<vector<int>> adj(n, vector<int>(n, 0));

    // Randomize node order for spanning tree generation
    vector<int> nodes(n);
    iota(nodes.begin(), nodes.end(), 0);

    for (int i = n - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        swap(nodes[i], nodes[j]);
    }

    // Build spanning tree to guarantee connectivity
    for (int i = 1; i < n; ++i) {
        int u = nodes[i];
        int target_index = rand() % i;
        int v = nodes[target_index];

        if (rand() % 2 == 0) {
            adj[u][v] = 1;
        } else {
            adj[v][u] = 1;
        }
    }

    // Add additional random edges based on density
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            bool edgeExists = (adj[i][j] == 1 || adj[j][i] == 1);

            if (!edgeExists) {
                double r = static_cast<double>(rand()) / RAND_MAX;
                if (r < density) {
                    if (rand() % 2 == 0) {
                        adj[i][j] = 1;
                    } else {
                        adj[j][i] = 1;
                    }
                }
            }
        }
    }

    return adj;
}

void GraphGenerator::saveGraphsToFile(const string& filename,
                                    const vector<vector<int>>& G,
                                    const vector<vector<int>>& H) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not create file " << filename << endl;
        return;
    }

    // Write first graph G (pattern)
    file << G.size() << endl;
    for (const auto& row : G) {
        for (size_t j = 0; j < row.size(); ++j) {
            file << row[j] << (j == row.size() - 1 ? "" : " ");
        }
        file << endl;
    }

    // Write second graph H (target)
    file << H.size() << endl;
    for (const auto& row : H) {
        for (size_t j = 0; j < row.size(); ++j) {
            file << row[j] << (j == row.size() - 1 ? "" : " ");
        }
        file << endl;
    }

    file.close();
}
