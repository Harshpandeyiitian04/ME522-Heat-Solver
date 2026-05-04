# Parallel 2D Steady-State Heat Conduction Solver
**ME-522: High-Performance Scientific Computing — IIT Mandi**  
**Student:** Harsh Sharma (B23441)

---

## Overview
A Fortran 90 implementation of a parallel 2D steady-state heat conduction solver
using the Finite Difference Method (FDM) and Jacobi iteration, parallelised with
OpenMP.

**Governing equation:** ∂²T/∂x² + ∂²T/∂y² = 0  
**Domain:** [0,1] × [0,1], uniform N×N interior grid  
**Boundary conditions:** T = 100°C (top), T = 0°C (bottom/left/right)  
**Method:** Jacobi iterative solver, convergence criterion ||ΔT||∞ < 1e-6

---

## Repository Structure

```
heat2d/
├── src/
│   ├── heat_serial.f90      # Serial solver + cache ordering test
│   ├── heat_parallel.f90    # OpenMP parallel solver
│   └── heat_verify.f90      # Verification against analytical solution
├── scripts/
│   ├── run_scaling.sh       # Automated strong-scaling study
│   └── plot_results.py      # All matplotlib plots
├── data/                    # CSV timing results + temperature field output
├── plots/                   # Generated figures (PDF)
├── report/
│   └── report.tex           # Full technical report (LaTeX)
└── Makefile                 # Build and run system
```

---

## Quick Start

### Prerequisites
```bash
# Ubuntu / Debian
sudo apt-get install gfortran python3-pip
pip3 install numpy matplotlib
```

### Build and Run
```bash
# Build everything
make all

# Run verification (confirms O(h^2) convergence)
make verify

# Run serial timing for N=128, 256, 512
./heat_serial 128
./heat_serial 256
./heat_serial 512

# Run parallel (e.g. 256×256 on 4 threads)
./heat_parallel 256 4

# Full automated scaling study
make run

# Run cache ordering comparison
./heat_serial 256 1

# Generate all plots
make plots

# Clean
make clean
```

---

## Makefile Targets

| Target | Description |
|--------|-------------|
| `make all` | Build serial + parallel + verify |
| `make serial` | Build serial solver |
| `make parallel` | Build OpenMP solver |
| `make verify` | Build and run verification |
| `make run` | Full scaling study (all grids, all thread counts) |
| `make cache` | Loop-ordering comparison |
| `make plots` | Generate all figures |
| `make clean` | Remove binaries and data files |

---

## Key Results

| N | Serial time | Speedup (8T) | Efficiency |
|---|-------------|--------------|------------|
| 128×128 | 1.04 s | 4.13× | 51.6% |
| 256×256 | 14.6 s | 4.13× | 51.6% |
| 512×512 | 65.8 s | 4.13× | 51.6% |

**Verification:** L₂ convergence ratio ≈ 3.92–4.79 (confirms O(h²) accuracy)

---

## Implementation Notes

- **Column-major loop order:** Inner loop iterates over `i` (first index),
  which is contiguous in Fortran's column-major layout — maximising cache line reuse.
- **OpenMP clauses:**  
  - `PRIVATE(i)` — avoids data race on inner loop counter  
  - `REDUCTION(max:err)` — combines per-thread maximum residuals  
  - `SCHEDULE(static)` — equal-sized row blocks for load balance
- **Two-array Jacobi:** Pure Jacobi uses separate `T` and `Tnew` arrays,
  making all N² updates data-independent and race-free.

---

## Compilation Flags

```bash
# Serial
gfortran -O2 -Wall -std=f2008 -o heat_serial src/heat_serial.f90

# Parallel
gfortran -O2 -Wall -std=f2008 -fopenmp -o heat_parallel src/heat_parallel.f90
```
