# Graph Algorithms

C++ implementation of graph subgraph isomorphism algorithms.

## Algorithms

- **Exact**: Finds optimal solution via exhaustive search (small graphs only)
- **Hungarian**: Polynomial approximation using assignment algorithm  
- **Greedy**: Fast approximation using degree matching

## Build

```bash
mkdir build && cd build
cmake ..
make
```

## Usage

Generate graphs:
```bash
./graph-generator 5 8    # G(5 nodes) vs H(8 nodes)
```

Run algorithms:
```bash
./graph-solver exact graphs.txt      # Optimal (slow)
./graph-solver hungarian graphs.txt  # Good approximation
```

## Test

```bash
./test.sh
```

Compares all three algorithms on small graphs. Exact works up to ~G(6), approximations scale to 1000+ nodes.