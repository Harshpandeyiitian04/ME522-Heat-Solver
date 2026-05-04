#!/bin/bash
#==============================================================================
# run_scaling.sh  –  Automated strong-scaling study
#
# Runs heat_serial and heat_parallel across all grid sizes and thread counts,
# appending one timing line per run to:
#   data/serial_timing.csv      (N, threads, wall_time, iterations)
#   data/parallel_timing.csv    (N, threads, wall_time, iterations)
#
# Header lines are written once at the start.
#==============================================================================

GRIDS="128 256 512"
THREADS="1 2 4 8"
DATADIR="data"

mkdir -p "$DATADIR"

# Write CSV headers (overwrite previous runs)
echo "N,threads,wall_time_s,iterations" > "$DATADIR/serial_timing.csv"
echo "N,threads,wall_time_s,iterations" > "$DATADIR/parallel_timing.csv"

echo "========================================================="
echo "  Strong Scaling Study"
echo "  Grids   : $GRIDS"
echo "  Threads : $THREADS"
echo "========================================================="

# ---- Serial baseline (1 thread) ------------------------------------------
for N in $GRIDS; do
    echo ""
    echo ">>> Serial  N=$N"
    ./heat_serial "$N"
done

# ---- Parallel (vary threads) ----------------------------------------------
for N in $GRIDS; do
    for T in $THREADS; do
        echo ""
        echo ">>> Parallel  N=$N  threads=$T"
        ./heat_parallel "$N" "$T"
    done
done

echo ""
echo "========================================================="
echo "  Scaling study complete."
echo "  Results in $DATADIR/serial_timing.csv"
echo "             $DATADIR/parallel_timing.csv"
echo "========================================================="
