#!/bin/bash

mkdir -p build
cd build
cmake ..
make -j4

if [ $? -ne 0 ]; then
    echo "Build failed"
    exit 1
fi

echo "Testing Graph Algorithms"
echo

for case in "3 5" "4 6" "5 8" "6 10"; do
    set -- $case
    echo "G($1) vs H($2):"

    ./graph-generator $1 $2 >/dev/null

    echo -n "  Exact: "
    exact_result=$(timeout 30s ./graph-solver exact graphs.txt 2>/dev/null)
    if [ $? -eq 124 ]; then
        echo "TIMEOUT"
    else
        exact_time=$(echo "$exact_result" | grep "Execution time:" | awk '{print $3}')
        exact_cost=$(echo "$exact_result" | grep "Minimal Extension Cost:" | awk '{print $4}')
        echo "${exact_time}ms (cost $exact_cost)"
    fi

    echo -n "  Hungarian: "
    result=$(./graph-solver hungarian graphs.txt 2>/dev/null)
    h_time=$(echo "$result" | grep "Hungarian time:" | awk '{for(i=1;i<=NF;i++) if($i~/^[0-9]/) print $i; exit}')
    h_edges=$(echo "$result" | grep "extension size" | awk '{print $5}')
    echo "${h_time}ms ($h_edges edges)"

    echo -n "  Greedy: "
    g_time=$(echo "$result" | grep "Greedy time:" | awk '{for(i=1;i<=NF;i++) if($i~/^[0-9]/) print $i; exit}')
    g_edges=$(echo "$result" | grep "Greedy edges added:" | awk '{print $4}')
    echo "${g_time}ms ($g_edges edges)"

    echo
done

echo "Done."
