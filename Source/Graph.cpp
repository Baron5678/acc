#include "Graph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <climits>
#include <functional>
#include <iostream>
#include <iomanip>
#include "HungarianAlgorithm.h"

#ifdef _WIN32
#include <windows.h>
#endif

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

void Graph::print() const {

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout << adj[i][j] << " ";
        }
        cout << endl;
    }
}

int Graph::edgeCount() const {
    int cnt = 0;
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            if (adj[i][j] != 0) cnt++;
    return cnt;
}


void Graph::printHighlighted(const Graph& other) const {
    using std::cout;
    using std::endl;

    const int nExt = this->size;
    const int nOrg = other.size;

    const int n = (nExt > nOrg) ? nExt : nOrg;

    const int W = 2;        
    const int GAP = 6;       
    const int LABELW = 6;    

#ifdef _WIN32
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

    CONSOLE_SCREEN_BUFFER_INFO csbi;
    WORD savedAttr = 0;
    if (hConsole != INVALID_HANDLE_VALUE && GetConsoleScreenBufferInfo(hConsole, &csbi)) {
        savedAttr = csbi.wAttributes;
    }
    else {
        hConsole = INVALID_HANDLE_VALUE;
    }

    const WORD GREEN = FOREGROUND_GREEN | FOREGROUND_INTENSITY;

    auto setColor = [&](WORD attr) {
        if (hConsole != INVALID_HANDLE_VALUE) SetConsoleTextAttribute(hConsole, attr);
        };
#else
    const char* GREEN = "\033[1;32m";
    const char* RESET = "\033[0m";
#endif

    auto getCell = [&](const Graph& g, int i, int j) -> int {
        if (i < 0 || j < 0) return 0;
        if (i >= g.size || j >= g.size) return 0;
        return g.adj[i][j];
        };

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << std::setw(W) << getCell(other, i, j);
        }

        cout << std::string(GAP, ' ');


        for (int j = 0; j < n; ++j) {
            const int org = getCell(other, i, j);
            const int ext = getCell(*this, i, j);

            const bool isNewEdge = (ext != 0 && org == 0);

#ifdef _WIN32
            if (isNewEdge) setColor(GREEN);
            cout << std::setw(W) << ext;
            if (isNewEdge) setColor(savedAttr);
#else
            if (isNewEdge) cout << GREEN;
            cout << std::setw(W) << ext;
            if (isNewEdge) cout << RESET;
#endif
        }

        cout << "\n";
}

#ifdef _WIN32
    if (hConsole != INVALID_HANDLE_VALUE) {
        SetConsoleTextAttribute(hConsole, savedAttr);
    }
#endif
}


int Graph::ComputeDistance(const Graph& other, const vector<int>& mapping) const {
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

pair<vector<int>, int> Graph::FindBestMapping(const Graph& target) const {
    vector<int> mapping(size, -1);
    vector<bool> usedH(target.size, false);
    
    function<pair<vector<int>, int>(int)> solve = [&](int vG) -> pair<vector<int>, int> {
        if (vG == this->size) {
            int distance = this->ComputeDistance(target, mapping);
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
