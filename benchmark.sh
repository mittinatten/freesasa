#!/bin/bash
# FreeSASA Parallel Scaling Benchmark
# Compares performance across thread counts for both algorithms

FREESASA="./freesasa"
PDB="../tests/data/1ubq.pdb"
REPEATS=5  # average over N runs

echo "================================================================"
echo "  FreeSASA Parallelized — Scaling Benchmark"
echo "================================================================"
echo ""
echo "PDB: $PDB"
echo "Repeats per measurement: $REPEATS"
echo ""

# Get system info
echo "System: $(nproc) cores available"
echo ""

# Function to benchmark
benchmark() {
    local alg=$1
    local threads=$2
    local alg_flag=""
    
    if [ "$alg" = "SR" ]; then
        alg_flag="--shrake-rupley --resolution=500"
    else
        alg_flag="--lee-richards --resolution=100"
    fi

    local total=0
    for r in $(seq 1 $REPEATS); do
        local t=$( { time $FREESASA $alg_flag --n-threads=$threads $PDB > /dev/null 2>&1 ; } 2>&1 | grep real | awk '{print $2}' | sed 's/m/*60+/' | sed 's/s//' | bc -l )
        total=$(echo "$total + $t" | bc -l)
    done
    local avg=$(echo "scale=4; $total / $REPEATS" | bc -l)
    echo "$avg"
}

# High-resolution benchmark for meaningful timing
echo "================================================================"
echo "  Shrake & Rupley (500 test points) — Higher workload"
echo "================================================================"
printf "%-10s %-12s %-10s\n" "Threads" "Time(s)" "Speedup"
printf "%-10s %-12s %-10s\n" "-------" "--------" "-------"

baseline=""
for t in 1 2 4 8 16; do
    time_val=$(benchmark "SR" $t)
    if [ "$t" -eq 1 ]; then
        baseline=$time_val
    fi
    speedup=$(echo "scale=2; $baseline / $time_val" | bc -l 2>/dev/null || echo "N/A")
    printf "%-10s %-12s %-10s\n" "$t" "${time_val}s" "${speedup}x"
done

echo ""
echo "================================================================"
echo "  Lee & Richards (100 slices) — Higher workload"
echo "================================================================"
printf "%-10s %-12s %-10s\n" "Threads" "Time(s)" "Speedup"
printf "%-10s %-12s %-10s\n" "-------" "--------" "-------"

baseline=""
for t in 1 2 4 8 16; do
    time_val=$(benchmark "LR" $t)
    if [ "$t" -eq 1 ]; then
        baseline=$time_val
    fi
    speedup=$(echo "scale=2; $baseline / $time_val" | bc -l 2>/dev/null || echo "N/A")
    printf "%-10s %-12s %-10s\n" "$t" "${time_val}s" "${speedup}x"
done

echo ""
echo "Benchmark complete."
