@echo off
REM ============================================================
REM  build.bat — Compile and run DNA Sequence Similarity
REM              Interactive Verification Mode
REM ============================================================

echo [1/2] Compiling with optimisations (-O2)...
g++ -O2 -std=c++17 -Wall -o dna_sim.exe main.cpp

IF ERRORLEVEL 1 (
    echo.
    echo ERROR: Compilation failed. Check that g++ is in your PATH.
    pause
    exit /b 1
)

echo [2/2] Build successful: dna_sim.exe
echo.
echo Launching Interactive Verification Mode...
echo.
dna_sim.exe --interactive

pause

