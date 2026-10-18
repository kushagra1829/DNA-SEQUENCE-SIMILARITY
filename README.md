# DNA Sequence Similarity — Performance Programming Project

## Project Structure

```
DNASequenceSimilarity/
├── main.cpp              ← All C++ code (algorithm + experiments)
├── build.bat             ← One-click build and run script (Windows)
├── graphs/
│   └── plot_results.py   ← Generates all 7 analysis graphs from CSV
├── results_benchmark.csv ← Generated after running dna_sim.exe
├── results_mutation.csv  ← Generated after running dna_sim.exe
└── calibration_log.csv   ← Generated after running dna_sim.exe
```

## Quick Start

### 1. Install MinGW-w64 (via MSYS2)
If you used the auto-installer, open a **new** PowerShell/terminal window after install, then install the compiler toolchain:
```
C:\msys64\msys2.exe
```
Inside the MSYS2 terminal run:
```bash
pacman -S --noconfirm mingw-w64-x86_64-gcc
```
Then add `C:\msys64\mingw64\bin` to your Windows PATH.

### 2. Compile
```powershell
cd DNASequenceSimilarity
g++ -O2 -std=c++17 -Wall -o dna_sim.exe main.cpp
```
Or simply double-click **build.bat**.

### 3. Run all experiments
```powershell
.\dna_sim.exe
```
This performs:
- **Calibration** — measures both k-mer strategies on a 50 000-bp test sequence, determines experimental threshold, writes `calibration_log.csv`
- **Correctness tests** — 9 test cases including hand-verified examples
- **Experiment A & B** — benchmarks lengths {1k, 10k, 100k, 1M} × k={1..32}, writes `results_benchmark.csv`
- **Experiment C** — mutation sensitivity at lengths 100 000 × k={3,8,16,24,32}, writes `results_mutation.csv`

### 4. Generate graphs
```powershell
pip install matplotlib pandas
python graphs\plot_results.py
```
Outputs to `graphs/`:
| File | Content |
|------|---------|
| `calibration_crossover.png` | Small vs large strategy timing per k |
| `exp_a_time_vs_length.png` | Time vs sequence length (log-log) |
| `exp_b_time_vs_k.png` | Time vs k for each sequence length |
| `phase_breakdown.png` | k-mer counting vs Jaccard time split |
| `exp_c_mutation_vs_similarity.png` | Similarity decay under mutation |
| `heatmap_time.png` | Colour heatmap of all [len × k] timings |
| `strategy_comparison.png` | Bar chart comparing two strategies |

## Algorithm Summary

### k-mer counting — Two Strategies

| Strategy | Data structure | Best for |
|----------|---------------|---------|
| `countKmers_small` | `std::unordered_map<string,int>` with pre-reserved buckets | Small k — few distinct keys, short strings, high frequency → O(1) hash lookup amortised |
| `countKmers_large` | `std::map<string,int>` (sorted tree) | Large k — many long unique keys → tree comparison short-circuits on prefix mismatch |

The **threshold** between the two strategies is determined **experimentally at startup** by the `calibrate()` function: it benchmarks both strategies across k=1..32 on a 50 000-bp sequence and selects the first k where the large strategy is ≥5% faster.

### Generalized Jaccard Index

```
J(A, B) = Σ min(c_A(x), c_B(x)) / Σ max(c_A(x), c_B(x))
```

- Score **1.0** → identical k-mer frequency distributions
- Score **0.0** → no k-mers in common
- Score **0..1** → partial similarity

### Complexity Estimates
- Building frequency map: O(n·k) character operations (each of the n-k+1 k-mers costs O(k) to hash/compare)
- Jaccard computation: O(d_A + d_B) hash-map operations, d = distinct k-mers
- Memory: O(d) entries × O(k) bytes per key

## Reproducibility
- All random sequences use a fixed seed (`std::mt19937` seeded with constants)
- Compiler: `g++ -O2 -std=c++17`
- Timing: `std::chrono::steady_clock` (monotonic, nanosecond resolution)
- Warmup run before each benchmark series to avoid cold-cache bias
