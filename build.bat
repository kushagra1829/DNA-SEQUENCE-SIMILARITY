@echo off
REM ============================================================
REM  build.bat — Compile and run the DNA Sequence Similarity
REM              project using MinGW-w64 g++
REM
REM  Prerequisites:
REM    Add C:\msys64\mingw64\bin to your PATH, or run this
REM    script from the MSYS2 MinGW 64-bit shell.
REM ============================================================

echo [1/3] Compiling with optimisations (-O2)...
g++ -O2 -std=c++17 -Wall -o dna_sim.exe main.cpp

IF ERRORLEVEL 1 (
    echo.
    echo ERROR: Compilation failed.  Check that g++ is in your PATH.
    echo   Typical path: C:\msys64\mingw64\bin\g++.exe
    pause
    exit /b 1
)

echo [2/3] Build successful.  Binary: dna_sim.exe
echo.
echo [3/3] Running experiments (this may take several minutes for len=1,000,000)...
dna_sim.exe

echo.
echo Done. CSV files are in the current directory.
echo Run  python graphs\plot_results.py  to generate graphs.
pause
