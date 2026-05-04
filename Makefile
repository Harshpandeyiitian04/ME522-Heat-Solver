#==============================================================================
# Makefile  –  ME-522 Project: Parallel 2D Heat Conduction Solver
# Author  : Harsh Sharma  (B23441)
# Course  : ME-522, IIT Mandi
#
# Targets
#   make serial    – build serial solver
#   make parallel  – build OpenMP parallel solver
#   make verify    – build and run verification against analytical solution
#   make all       – build serial + parallel + verify
#   make run       – full scaling study (calls run_scaling.sh)
#   make cache     – run loop-ordering comparison
#   make plots     – generate all Python plots
#   make clean     – remove binaries and output data
#==============================================================================

FC       = gfortran
FFLAGS   = -O2 -Wall -std=f2008
OMPFLAGS = -fopenmp
SRCDIR   = src
DATADIR  = data
PLOTDIR  = plots

# Grid sizes and thread counts used in the scaling study
GRIDS    = 128 256 512
THREADS  = 1 2 4 8

.PHONY: all serial parallel verify run cache plots clean help

# ---- Default target -------------------------------------------------------
all: serial parallel verify

# ---- Build targets --------------------------------------------------------
serial: $(SRCDIR)/heat_serial.f90
	@mkdir -p $(DATADIR)
	$(FC) $(FFLAGS) -o heat_serial $<
	@echo ">>> heat_serial built"

parallel: $(SRCDIR)/heat_parallel.f90
	@mkdir -p $(DATADIR)
	$(FC) $(FFLAGS) $(OMPFLAGS) -o heat_parallel $<
	@echo ">>> heat_parallel built"

verify: $(SRCDIR)/heat_verify.f90
	@mkdir -p $(DATADIR)
	$(FC) $(FFLAGS) -o heat_verify $<
	@rm -f $(DATADIR)/verify_errors.csv
	./heat_verify
	@echo ">>> Verification complete"

# ---- Run targets ----------------------------------------------------------
run: serial parallel
	@bash scripts/run_scaling.sh

cache: serial
	@echo ">>> Running loop-ordering (cache) comparison..."
	@rm -f $(DATADIR)/cache_timing.csv
	@for N in 256 512; do \
	    echo "  N=$$N ..."; \
	    ./heat_serial $$N 1; \
	done

plots: scripts/plot_results.py
	@mkdir -p $(PLOTDIR)
	python3 scripts/plot_results.py
	@echo ">>> Plots written to $(PLOTDIR)/"

# ---- Utility targets ------------------------------------------------------
clean:
	rm -f heat_serial heat_parallel heat_verify *.mod
	rm -f $(DATADIR)/serial_timing.csv $(DATADIR)/parallel_timing.csv
	rm -f $(DATADIR)/verify_errors.csv $(DATADIR)/cache_timing.csv
	rm -f $(DATADIR)/temp_*.dat

help:
	@echo ""
	@echo "  make all       build serial + parallel + verify"
	@echo "  make run       run full scaling experiments"
	@echo "  make cache     run loop-ordering comparison"
	@echo "  make plots     generate publication plots (requires matplotlib)"
	@echo "  make clean     remove binaries and data files"
	@echo ""
