#include "GraphGenerator.h"
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main(int argc, char* argv[]) {
    srand(static_cast<unsigned int>(time(0)));

    int N_G = 5;
    int N_H = 8;

    if (argc == 3) {
        try {
            N_G = stoi(argv[1]);
            N_H = stoi(argv[2]);
        }
        catch (...) {
            cerr << "Error: Invalid arguments. Usage: " << argv[0] << " {size V(G)} {size V(H)}" << endl;
            return 1;
        }
    }
    else if (argc > 1) {
        cout << "Usage: " << argv[0] << " {size V(G)} {size V(H)}" << endl;
        cout << "Using default sizes: " << N_G << ", " << N_H << endl;
    }

    if (N_H < N_G) {
        cout << "Warning: Target H (" << N_H << ") is smaller than Pattern G (" << N_G << ")." << endl;
    }

    double density_G = 0.8;
    double density_H = 0.3;

    auto G = GraphGenerator::generateConnectedGraph(N_G, density_G);
    auto H = GraphGenerator::generateConnectedGraph(N_H, density_H);

    string filename = "graphs.txt";
    GraphGenerator::saveGraphsToFile(filename, G, H);

    cout << "Generated graphs: G(" << N_G << ") and H(" << N_H << ")" << endl;
    cout << "Saved to: " << filename << endl;

    return 0;
}
