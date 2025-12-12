#include "../include/HungarianAlgorithm.h"
#include <climits>
#include <vector>

using namespace std;

HungarianAlgorithm::HungarianAlgorithm(int size) : n(size) {
    cost_matrix.assign(n, vector<int>(n, 0));
}

void HungarianAlgorithm::setCost(int i, int j, int cost) {
    cost_matrix[i][j] = cost;
}

vector<int> HungarianAlgorithm::findMinCostAssignment() {
    return solve();
}

vector<int> HungarianAlgorithm::solve() {
    // Kuhn-Munkres algorithm for minimum cost bipartite matching
    vector<vector<int>> matrix = cost_matrix;
    vector<int> u(n + 1), v(n + 1), p(n + 1), way(n + 1);

    for (int i = 1; i <= n; ++i) {
        // Process each row to build optimal assignment
        p[0] = i;
        int j0 = 0;
        vector<int> minv(n + 1, INT_MAX);
        vector<bool> used(n + 1, false);

        do {
            // Find augmenting path using modified Dijkstra
            used[j0] = true;
            int i0 = p[j0], delta = INT_MAX, j1;

            // Update reduced costs for all unmatched columns
            for (int j = 1; j <= n; ++j) {
                if (!used[j]) {
                    int cur = matrix[i0 - 1][j - 1] - u[i0] - v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }

            // Update dual variables (potentials)
            for (int j = 0; j <= n; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }

            j0 = j1;
        } while (p[j0] != 0);

        // Update assignment along augmenting path
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }

    // Convert internal representation to result format
    vector<int> result(n);
    for (int j = 1; j <= n; ++j) {
        if (p[j] != 0) {
            result[p[j] - 1] = j - 1;
        }
    }
    return result;
}
