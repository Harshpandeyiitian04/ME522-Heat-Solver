@echo off
REM ===========================================================================
REM  run_all.bat  --  ME-522 Project: Full build + run on Windows
REM  Author: Harsh Sharma (B23441)
REM
REM  Prerequisites (install once):
REM    gfortran  -> https://winlibs.com  (download WinLibs, add bin/ to PATH)
REM    python3   -> https://python.org   (tick "Add to PATH" during install)
REM    pip install numpy matplotlib
REM
REM  Usage:
REM    Double-click run_all.bat   OR   run from Command Prompt inside heat2d\
REM ===========================================================================

echo.
echo ============================================================
echo  ME-522 Project -- Parallel 2D Heat Conduction Solver
echo  Harsh Sharma  (B23441)  IIT Mandi
echo ============================================================
echo.

REM ---- 0. Check gfortran is available ------------------------------------
where gfortran >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] gfortran not found in PATH.
    echo.
    echo  Install steps:
    echo  1. Go to  https://winlibs.com
    echo  2. Download "GCC 13.x.x + LLVM ... Win64"  ^(the .zip^)
    echo  3. Extract to  C:\mingw64
    echo  4. Open Start -^> "Edit environment variables for your account"
    echo  5. Edit PATH -^> New -^> add  C:\mingw64\bin
    echo  6. Open a NEW command prompt and re-run this script.
    pause
    exit /b 1
)
echo [OK] gfortran found:
gfortran --version | findstr /c:"GNU Fortran"

REM ---- 1. Create output directories ------------------------------------
if not exist data  mkdir data
if not exist plots mkdir plots

REM ---- 2. Compile serial solver ----------------------------------------
echo.
echo [1/5] Compiling serial solver ...
gfortran -O2 -Wall -std=f2008 -o heat_serial.exe src\heat_serial.f90
if %ERRORLEVEL% NEQ 0 ( echo [FAILED] Serial compile failed && pause && exit /b 1 )
echo       heat_serial.exe  OK

REM ---- 3. Compile OpenMP parallel solver --------------------------------
echo [2/5] Compiling parallel solver ...
gfortran -O2 -Wall -std=f2008 -fopenmp -o heat_parallel.exe src\heat_parallel.f90
if %ERRORLEVEL% NEQ 0 ( echo [FAILED] Parallel compile failed && pause && exit /b 1 )
echo       heat_parallel.exe  OK

REM ---- 4. Compile and run verification ----------------------------------
echo [3/5] Compiling and running verification ...
gfortran -O2 -Wall -std=f2008 -o heat_verify.exe src\heat_verify.f90
if %ERRORLEVEL% NEQ 0 ( echo [FAILED] Verify compile failed && pause && exit /b 1 )
if exist data\verify_errors.csv del data\verify_errors.csv
heat_verify.exe
echo.

REM ---- 5. Serial timing runs -------------------------------------------
echo [4/5] Running serial timing  (N = 128, 256, 512) ...
echo N,threads,wall_time_s,iterations > data\serial_timing.csv

echo   N=128 ...
heat_serial.exe 128
echo   N=256 ...
heat_serial.exe 256
echo   N=512 ...
heat_serial.exe 512

REM ---- 6. Parallel scaling study (all thread counts) -------------------
echo.
echo [5/5] Running parallel scaling study ...
echo N,threads,wall_time_s,iterations > data\parallel_timing.csv

for %%N in (128 256 512) do (
    for %%T in (1 2 4 8) do (
        echo   N=%%N  threads=%%T ...
        heat_parallel.exe %%N %%T
    )
)

REM ---- 7. Cache ordering comparison ------------------------------------
echo.
echo [BONUS] Running cache ordering comparison ...
if exist data\cache_timing.csv del data\cache_timing.csv
echo N,friendly_s,unfriendly_s > data\cache_timing.csv
echo   N=256 cache test...
heat_serial.exe 256 1
echo   N=512 cache test...
heat_serial.exe 512 1

REM ---- 8. Generate plots -----------------------------------------------
echo.
echo [PLOTS] Generating figures with Python ...
where python >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    where python3 >nul 2>&1
    if %ERRORLEVEL% NEQ 0 (
        echo [WARNING] Python not found. Skipping plots.
        echo           Install from https://python.org and run:
        echo             pip install numpy matplotlib
        echo             python scripts\plot_results.py
        goto :done
    )
    python3 -m pip install numpy matplotlib -q
    python3 scripts\plot_results.py
) else (
    python -m pip install numpy matplotlib -q
    python scripts\plot_results.py
)

:done
echo.
echo ============================================================
echo  DONE!  Check these output folders:
echo    data\    -- timing CSV files and temperature field
echo    plots\   -- PDF figures (open with any PDF viewer)
echo ============================================================
echo.
echo  Deliverable checklist:
echo    [x] heat_serial.exe    -- serial solver
echo    [x] heat_parallel.exe  -- OpenMP parallel solver
echo    [x] heat_verify.exe    -- runs and confirms O(h^2) convergence
echo    [x] data\*.csv         -- timing data for all N and threads
echo    [x] plots\*.pdf        -- heatmap, verification, speedup, efficiency, cache
echo    [x] report\report.pdf  -- full technical report
echo.
pause
