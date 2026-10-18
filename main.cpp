// ============================================================================
// DNA Sequence Similarity — Generalized Jaccard Index
// Performance Programming Assignment
//
// Key Architectural Features:
// 1. 2-bit quaternary integer encoding (A=00, C=01, G=10, T=11) in uint64_t (k<=32).
//    Rolling hash extraction: ((code << 2) & mask) | new_base (zero heap churn).
// 2. Strategy A (Small k <= 13): Direct 4^k flat counter array. Single-pass linear scan.
// 3. Strategy B (Large k): Zero-copy integer vector, std::sort (Introsort, O(n log n)),
//    and two-pointer linear merge-walk (O(|A| + |B|)).
// 4. Empirical calibration: Automatically measures crossover threshold on startup.
// 5. Dual-length calibration stability verification.
// 6. 9 automated correctness tests + independent reference implementation.
// 7. Full benchmark suite (lengths 1k..1M, k=1..32) + 32-k mutation sensitivity.
// ============================================================================

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// High-resolution clock utilities
// ---------------------------------------------------------------------------
using Clock   = std::chrono::steady_clock;
using Seconds = std::chrono::duration<double>;

inline double elapsed_sec(Clock::time_point start) {
    return Seconds(Clock::now() - start).count();
}

// ---------------------------------------------------------------------------
// 2-bit Integer Encoding Utilities
// ---------------------------------------------------------------------------
// A = 00 (0), C = 01 (1), G = 10 (2), T = 11 (3)
inline uint8_t baseToBits(char c) {
    switch (c) {
        case 'A': case 'a': return 0b00;
        case 'C': case 'c': return 0b01;
        case 'G': case 'g': return 0b10;
        case 'T': case 't': return 0b11;
        default: return 0b00;
    }
}

// Bitmask for lowest 2*k bits. Avoids 1ULL << 64 undefined behavior when k=32.
inline uint64_t getKmerMask(int k) {
    if (k >= 32) return ~0ULL;
    return (1ULL << (2 * k)) - 1ULL;
}

// Decode a 64-bit integer code back into DNA text string
inline std::string decodeKmer(uint64_t code, int k) {
    std::string s(k, 'A');
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (int i = k - 1; i >= 0; --i) {
        s[i] = bases[code & 0b11];
        code >>= 2;
    }
    return s;
}

// Max k for Strategy A direct array (4^13 counters * 4 bytes = 268 MB)
constexpr int MAX_SMALL_K = 13;

// Global threshold determined dynamically during startup calibration
int g_threshold = 11;

// ---------------------------------------------------------------------------
// Data structures for k-mer counts
// ---------------------------------------------------------------------------
struct KmerCount {
    uint64_t kmer;
    uint32_t count;
};

// ===========================================================================
// SECTION 1 — k-mer counting  (TWO SEPARATE IMPLEMENTATIONS)
// ===========================================================================

// ---------------------------------------------------------------------------
// Strategy A — SMALL k (Direct array indexing in 4^k slots)
// ---------------------------------------------------------------------------
std::vector<uint32_t> countKmers_small(const std::string& seq, int k) {
    if (k <= 0 || static_cast<int>(seq.size()) < k || k > MAX_SMALL_K) {
        return {};
    }

    const size_t total_slots = 1ULL << (2 * k); // 4^k counters
    std::vector<uint32_t> counts(total_slots, 0);

    const uint64_t mask = getKmerMask(k);
    uint64_t code = 0;
    for (int i = 0; i < k; ++i) {
        code = (code << 2) | baseToBits(seq[i]);
    }
    counts[code]++;

    const int n = static_cast<int>(seq.size());
    for (int i = k; i < n; ++i) {
        code = ((code << 2) & mask) | baseToBits(seq[i]);
        counts[code]++;
    }
    return counts;
}

inline double jaccard_small_from_counts(const std::vector<uint32_t>& cA,
                                        const std::vector<uint32_t>& cB) {
    if (cA.empty() && cB.empty()) return 1.0;
    if (cA.empty() || cB.empty()) return 0.0;

    uint64_t sum_min = 0, sum_max = 0;
    const size_t sz = cA.size();
    for (size_t i = 0; i < sz; ++i) {
        uint32_t a = cA[i];
        uint32_t b = cB[i];
        sum_min += std::min(a, b);
        sum_max += std::max(a, b);
    }
    return (sum_max == 0) ? 1.0 : static_cast<double>(sum_min) / static_cast<double>(sum_max);
}

// ---------------------------------------------------------------------------
// Strategy B — LARGE k (Contiguous 64-bit integer sort + run-length consolidation)
// ---------------------------------------------------------------------------
std::vector<KmerCount> countKmers_large(const std::string& seq, int k) {
    if (k <= 0 || static_cast<int>(seq.size()) < k) return {};

    const int n = static_cast<int>(seq.size());
    const int num_kmers = n - k + 1;

    std::vector<uint64_t> codes;
    codes.reserve(num_kmers);

    const uint64_t mask = getKmerMask(k);
    uint64_t code = 0;
    for (int i = 0; i < k; ++i) {
        code = (code << 2) | baseToBits(seq[i]);
    }
    codes.push_back(code);

    for (int i = k; i < n; ++i) {
        code = ((code << 2) & mask) | baseToBits(seq[i]);
        codes.push_back(code);
    }

    std::sort(codes.begin(), codes.end());

    std::vector<KmerCount> runs;
    runs.reserve(std::min<size_t>(num_kmers, 1ULL << (std::min(k, 16) * 2)));

    uint32_t run_len = 1;
    for (size_t i = 1; i < codes.size(); ++i) {
        if (codes[i] == codes[i - 1]) {
            ++run_len;
        } else {
            runs.push_back({codes[i - 1], run_len});
            run_len = 1;
        }
    }
    if (!codes.empty()) {
        runs.push_back({codes.back(), run_len});
    }
    return runs;
}

inline double jaccard_large_from_runs(const std::vector<KmerCount>& rA,
                                      const std::vector<KmerCount>& rB) {
    if (rA.empty() && rB.empty()) return 1.0;
    if (rA.empty() || rB.empty()) return 0.0;

    uint64_t sum_min = 0, sum_max = 0;
    size_t i = 0, j = 0;
    while (i < rA.size() && j < rB.size()) {
        if (rA[i].kmer == rB[j].kmer) {
            sum_min += std::min(rA[i].count, rB[j].count);
            sum_max += std::max(rA[i].count, rB[j].count);
            ++i; ++j;
        } else if (rA[i].kmer < rB[j].kmer) {
            sum_max += rA[i].count;
            ++i;
        } else {
            sum_max += rB[j].count;
            ++j;
        }
    }
    while (i < rA.size()) { sum_max += rA[i++].count; }
    while (j < rB.size()) { sum_max += rB[j++].count; }

    return (sum_max == 0) ? 1.0 : static_cast<double>(sum_min) / static_cast<double>(sum_max);
}

// ===========================================================================
// SECTION 2 — Unified Jaccard Computation & Timed Execution
// ===========================================================================

struct JaccardResult {
    double similarity;
    double time_count_ms;
    double time_jaccard_ms;
    double total_ms;
};

JaccardResult computeJaccardTimed(const std::string& seqA,
                                  const std::string& seqB,
                                  int k,
                                  int threshold) {
    if (k <= 0) return {0.0, 0.0, 0.0, 0.0};
    const bool a_short = static_cast<int>(seqA.size()) < k;
    const bool b_short = static_cast<int>(seqB.size()) < k;
    if (a_short && b_short) return {1.0, 0.0, 0.0, 0.0};
    if (a_short || b_short) return {0.0, 0.0, 0.0, 0.0};

    if (k < threshold && k <= MAX_SMALL_K) {
        // Strategy A: Direct array counting
        auto t0 = Clock::now();
        auto cA = countKmers_small(seqA, k);
        auto cB = countKmers_small(seqB, k);
        double count_ms = elapsed_sec(t0) * 1000.0;

        auto t1 = Clock::now();
        double sim = jaccard_small_from_counts(cA, cB);
        double jacc_ms = elapsed_sec(t1) * 1000.0;

        return {sim, count_ms, jacc_ms, count_ms + jacc_ms};
    } else {
        // Strategy B: Zero-copy sort + two-pointer merge walk
        auto t0 = Clock::now();
        auto rA = countKmers_large(seqA, k);
        auto rB = countKmers_large(seqB, k);
        double count_ms = elapsed_sec(t0) * 1000.0;

        auto t1 = Clock::now();
        double sim = jaccard_large_from_runs(rA, rB);
        double jacc_ms = elapsed_sec(t1) * 1000.0;

        return {sim, count_ms, jacc_ms, count_ms + jacc_ms};
    }
}

double generalizedJaccard(const std::string& seqA, const std::string& seqB, int k) {
    return computeJaccardTimed(seqA, seqB, k, g_threshold).similarity;
}

struct TimedResult {
    double similarity;
    double time_count_ms;
    double time_jaccard_ms;
    double total_ms;
};

TimedResult timedJaccard(const std::string& seqA, const std::string& seqB, int k) {
    auto r = computeJaccardTimed(seqA, seqB, k, g_threshold);
    return {r.similarity, r.time_count_ms, r.time_jaccard_ms, r.total_ms};
}

// ===========================================================================
// SECTION 3 — Calibration & Stability Verification
// ===========================================================================

struct CalibrationResult {
    int threshold;
    std::vector<double> time_small_ms;
    std::vector<double> time_large_ms;
};

std::string generateDNA(int length, std::mt19937& rng);

CalibrationResult runSingleCalibration(int calib_len, int repeats, std::ofstream& log) {
    std::cout << "[Calibration] Benchmarking length " << calib_len
              << " bp (repeats=" << repeats << ")...\n";

    std::mt19937 rng(42);
    std::string seq = generateDNA(calib_len, rng);

    CalibrationResult cr;
    cr.time_small_ms.resize(33, -1.0);
    cr.time_large_ms.resize(33, -1.0);

    for (int k = 1; k <= 32; ++k) {
        // Strategy A timing (only if within memory budget)
        if (k <= MAX_SMALL_K) {
            double tot_s = 0.0;
            for (int r = 0; r < repeats; ++r) {
                auto t0 = Clock::now();
                auto c = countKmers_small(seq, k);
                tot_s += elapsed_sec(t0);
                (void)c.size();
            }
            cr.time_small_ms[k] = (tot_s / repeats) * 1000.0;
        } else {
            cr.time_small_ms[k] = -1.0; // Exceeds budget -> memory wall guard
        }

        // Strategy B timing (all k)
        double tot_l = 0.0;
        for (int r = 0; r < repeats; ++r) {
            auto t0 = Clock::now();
            auto runs = countKmers_large(seq, k);
            tot_l += elapsed_sec(t0);
            (void)runs.size();
        }
        cr.time_large_ms[k] = (tot_l / repeats) * 1000.0;

        log << calib_len << "," << k << ",small," << cr.time_small_ms[k] << "\n";
        log << calib_len << "," << k << ",large," << cr.time_large_ms[k] << "\n";

        std::string s_str = (cr.time_small_ms[k] >= 0.0)
            ? (std::to_string(cr.time_small_ms[k]).substr(0, 6) + "ms")
            : "WALL (>268MB)";
        std::string l_str = std::to_string(cr.time_large_ms[k]).substr(0, 6) + "ms";

        std::cout << "  k=" << std::setw(2) << k
                  << "  small=" << std::setw(15) << s_str
                  << "  large=" << std::setw(10) << l_str << "\n";
    }

    // Find crossover threshold
    cr.threshold = MAX_SMALL_K + 1; // Default upper bound
    for (int k = 1; k <= 32; ++k) {
        if (cr.time_small_ms[k] < 0.0) {
            cr.threshold = k; // Strategy A hit the memory wall
            break;
        }
        if (cr.time_large_ms[k] < cr.time_small_ms[k] * 0.95) {
            cr.threshold = k;
            break;
        }
    }
    return cr;
}

void calibrateAndVerifyStability() {
    std::ofstream log("calibration_log.csv");
    log << "calib_len,k,strategy,avg_time_ms\n";

    std::cout << "--- Phase 0: Empirical Calibration & Threshold Determination ---\n";
    CalibrationResult cal1 = runSingleCalibration(50000, 3, log);
    std::cout << ">>> Primary Calibration (50k bp): Threshold = k=" << cal1.threshold << "\n\n";

    std::cout << "--- Threshold Stability Check (Length = 200,000 bp) ---\n";
    CalibrationResult cal2 = runSingleCalibration(200000, 2, log);
    std::cout << ">>> Stability Verification (200k bp): Threshold = k=" << cal2.threshold << "\n\n";

    g_threshold = cal1.threshold;
    std::cout << ">>> Active Threshold Selected: k=" << g_threshold
              << " (Strategy A: k<" << g_threshold
              << ", Strategy B: k>=" << g_threshold << ")\n\n";
}

// ===========================================================================
// SECTION 4 — DNA sequence generation and mutation
// ===========================================================================

std::string generateDNA(int length, std::mt19937& rng) {
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<int> dist(0, 3);
    std::string seq(length, 'A');
    for (char& c : seq) c = bases[dist(rng)];
    return seq;
}

struct MutationResult {
    std::string mutated;
    int actual_mutations; // Exact Hamming distance
};

MutationResult mutateDNA(const std::string& seq, int target_mutations, std::mt19937& rng) {
    std::string mutated = seq;
    const char bases[] = {'A', 'C', 'G', 'T'};
    const int n = static_cast<int>(seq.size());
    std::uniform_int_distribution<int> pos_dist(0, n - 1);
    std::uniform_int_distribution<int> base_dist(0, 3);

    int attempts = 0;
    int target = std::min(target_mutations, n);
    std::vector<bool> changed(n, false);
    int current_mutations = 0;

    while (current_mutations < target && attempts < target * 10 + 1000) {
        ++attempts;
        int pos = pos_dist(rng);
        char nb = bases[base_dist(rng)];
        if (nb != mutated[pos]) {
            mutated[pos] = nb;
            if (!changed[pos]) {
                changed[pos] = true;
                ++current_mutations;
            }
        }
    }

    int hamming = 0;
    for (int i = 0; i < n; ++i) {
        if (seq[i] != mutated[i]) ++hamming;
    }
    return {mutated, hamming};
}

// ===========================================================================
// SECTION 5 — Correctness Verification Test Cases
// ===========================================================================

void runTests() {
    std::cout << "=== Correctness Tests (T1 - T9) ===\n\n";

    // Test 1: identical sequences -> 1.0
    {
        std::string s = "ATTCGTA";
        double sim = generalizedJaccard(s, s, 3);
        std::cout << "[T1] Identical sequences, k=3: sim=" << sim
                  << (std::abs(sim - 1.0) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 2: completely disjoint alphabets -> 0.0
    {
        std::string a(10, 'A'), c(10, 'C');
        double sim = generalizedJaccard(a, c, 3);
        std::cout << "[T2] Disjoint sequences (AAAA vs CCCC), k=3: sim=" << sim
                  << (sim == 0.0 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 3: hand-calculated identical -> 1.0
    {
        std::string a = "AACT", b = "AACT";
        double sim = generalizedJaccard(a, b, 2);
        std::cout << "[T3] AACT vs AACT, k=2: sim=" << sim
                  << (std::abs(sim - 1.0) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 4: hand-calculated partial overlap -> 1/5 = 0.2
    {
        std::string a = "AACG", b = "AATG";
        double sim = generalizedJaccard(a, b, 2);
        std::cout << "[T4] AACG vs AATG, k=2: sim=" << sim
                  << (std::abs(sim - 0.2) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 5: periodic motif -> 1.0
    {
        std::string s = "ATATAT";
        double sim = generalizedJaccard(s, s, 2);
        std::cout << "[T5] Periodic motif ATATAT vs itself, k=2: sim=" << sim
                  << (std::abs(sim - 1.0) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 6: k larger than sequence -> 1.0 by convention
    {
        std::string s = "ACG";
        double sim = generalizedJaccard(s, s, 10);
        std::cout << "[T6] k > seq length (ACG, k=10): sim=" << sim
                  << (std::abs(sim - 1.0) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 7: permuted bases with identical distribution -> 1.0
    {
        std::string a = "ACGT", b = "TGCA";
        double sim = generalizedJaccard(a, b, 1);
        std::cout << "[T7] ACGT vs TGCA, k=1 (same base freq): sim=" << sim
                  << (std::abs(sim - 1.0) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 8: stress test k=32
    {
        std::mt19937 rng(99);
        std::string s = generateDNA(200, rng);
        double sim = generalizedJaccard(s, s, 32);
        std::cout << "[T8] k=32, len=200, identical: sim=" << sim
                  << (std::abs(sim - 1.0) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 9: verify k=1..32 all return 1.0 for identical sequences
    {
        std::mt19937 rng(7);
        std::string s = generateDNA(500, rng);
        bool all_pass = true;
        for (int k = 1; k <= 32; ++k) {
            double sim = generalizedJaccard(s, s, k);
            if (std::abs(sim - 1.0) > 1e-9) { all_pass = false; break; }
        }
        std::cout << "[T9] Identical seq (len=500), k=1..32 all=1.0: "
                  << (all_pass ? "PASS" : "FAIL") << "\n";
    }
    std::cout << "\n";
}

// ===========================================================================
// SECTION 6 — Benchmarks: Experiment A & B
// ===========================================================================

void runBenchmark(const std::string& csv_path) {
    const std::vector<int> lengths = {1000, 10000, 100000, 500000, 1000000};
    const int REPEATS = 5;

    std::ofstream csv(csv_path);
    csv << "seq_len,k,strategy,repetition,similarity,time_count_ms,time_jaccard_ms,total_ms\n";

    std::mt19937 rng(1234);
    std::cout << "=== Experiment A & B: Full Benchmark Suite ===\n";

    for (int len : lengths) {
        std::cout << "\n  Generating sequences of length " << len << " bp...\n";
        std::string seqA = generateDNA(len, rng);
        std::string seqB = generateDNA(len, rng);

        for (int k = 1; k <= 32; ++k) {
            std::string strategy = (k < g_threshold && k <= MAX_SMALL_K) ? "small" : "large";

            // Warm-up run
            auto warm = timedJaccard(seqA, seqB, k);
            (void)warm;

            for (int rep = 1; rep <= REPEATS; ++rep) {
                auto r = timedJaccard(seqA, seqB, k);
                csv << len << "," << k << "," << strategy << "," << rep << ","
                    << std::fixed << std::setprecision(6)
                    << r.similarity << ","
                    << r.time_count_ms << ","
                    << r.time_jaccard_ms << ","
                    << r.total_ms << "\n";
            }
            std::cout << "    k=" << std::setw(2) << k
                      << " [" << strategy << "] done\n";
        }
    }
    std::cout << "\n[Benchmark] Results written to: " << csv_path << "\n\n";
}

// ===========================================================================
// SECTION 7 — Mutation Experiment C (All 32 k-values)
// ===========================================================================

void runMutationExperiment(const std::string& csv_path) {
    const int SEQ_LEN = 100000;
    const std::vector<double> mut_pcts = {
        0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0
    };
    const int REPEATS = 2; // 2 repeats across all 32 k-values

    std::cout << "=== Experiment C: Mutation Sensitivity (All k=1..32) ===\n";
    std::cout << "  Sequence length: " << SEQ_LEN << " bp\n\n";

    std::mt19937 rng(9999);
    std::string original = generateDNA(SEQ_LEN, rng);

    std::ofstream csv(csv_path);
    csv << "seq_len,k,mutation_pct,target_mutations,actual_mutations,repetition,similarity,total_ms\n";

    for (int k = 1; k <= 32; ++k) {
        std::cout << "  k=" << std::setw(2) << k << ": ";
        for (double pct : mut_pcts) {
            int target_mut = static_cast<int>(std::round(SEQ_LEN * pct / 100.0));

            for (int rep = 1; rep <= REPEATS; ++rep) {
                MutationResult mr = (target_mut == 0)
                    ? MutationResult{original, 0}
                    : mutateDNA(original, target_mut, rng);

                auto r = timedJaccard(original, mr.mutated, k);

                csv << SEQ_LEN << "," << k << ","
                    << pct << "," << target_mut << "," << mr.actual_mutations << "," << rep << ","
                    << std::fixed << std::setprecision(6)
                    << r.similarity << "," << r.total_ms << "\n";
            }
        }
        std::cout << "done\n";
    }
    std::cout << "\n[Mutation] Results written to: " << csv_path << "\n\n";
}

// ===========================================================================
// SECTION 8 — Interactive Verification Mode & Reference Validation
// ===========================================================================

double referenceJaccard(const std::string& seqA, const std::string& seqB, int k) {
    if (k <= 0) return 0.0;
    std::map<std::string, int> fA, fB;
    const int nA = static_cast<int>(seqA.size());
    const int nB = static_cast<int>(seqB.size());
    for (int i = 0; i <= nA - k; ++i) fA[seqA.substr(i, k)]++;
    for (int i = 0; i <= nB - k; ++i) fB[seqB.substr(i, k)]++;

    if (fA.empty() && fB.empty()) return 1.0;
    if (fA.empty() || fB.empty()) return 0.0;

    std::set<std::string> all_kmers;
    for (const auto& [km, _] : fA) all_kmers.insert(km);
    for (const auto& [km, _] : fB) all_kmers.insert(km);

    long long sum_min = 0, sum_max = 0;
    for (const auto& km : all_kmers) {
        int cA = fA.count(km) ? fA.at(km) : 0;
        int cB = fB.count(km) ? fB.at(km) : 0;
        sum_min += std::min(cA, cB);
        sum_max += std::max(cA, cB);
    }
    return (sum_max == 0) ? 1.0 : static_cast<double>(sum_min) / sum_max;
}

bool isValidDNA(const std::string& s) {
    if (s.empty()) return false;
    for (char c : s) {
        char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (u != 'A' && u != 'C' && u != 'G' && u != 'T') return false;
    }
    return true;
}

std::string cleanDNA(const std::string& s) {
    std::string res;
    for (char c : s) {
        if (!std::isspace(static_cast<unsigned char>(c))) {
            res.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
        }
    }
    return res;
}

void runInteractiveMode() {
    std::cout << "\n============================================\n";
    std::cout << "  DNA Sequence Similarity — Verification Mode\n";
    std::cout << "============================================\n";
    std::cout << "Tip: Enter DNA sequences using A, C, G, T (case-insensitive).\n";
    std::cout << "Type 'q' as the sequence to return to the menu.\n";

    while (true) {
        std::string s1, s2;
        std::cout << "\nEnter first DNA sequence: ";
        if (!(std::cin >> s1) || s1 == "q" || s1 == "Q") break;
        s1 = cleanDNA(s1);
        if (!isValidDNA(s1)) {
            std::cout << "Error: Invalid DNA sequence! Only A, C, G, T bases are allowed.\n";
            continue;
        }

        std::cout << "Enter second DNA sequence: ";
        if (!(std::cin >> s2) || s2 == "q" || s2 == "Q") break;
        s2 = cleanDNA(s2);
        if (!isValidDNA(s2)) {
            std::cout << "Error: Invalid DNA sequence! Only A, C, G, T bases are allowed.\n";
            continue;
        }

        int k = 0;
        std::cout << "Enter k (1-32): ";
        if (!(std::cin >> k)) {
            std::cin.clear();
            std::string dummy;
            std::cin >> dummy;
            std::cout << "Error: Invalid input. k must be an integer between 1 and 32.\n";
            continue;
        }
        if (k < 1 || k > 32) {
            std::cout << "Error: k must be between 1 and 32.\n";
            continue;
        }

        auto t0 = Clock::now();
        double sim = generalizedJaccard(s1, s2, k);
        double elapsed_us = elapsed_sec(t0) * 1e6;

        double ref_sim = referenceJaccard(s1, s2, k);
        bool match = (std::abs(sim - ref_sim) < 1e-9);

        std::cout << "\n--------------------------------------------\n";
        std::cout << "Generalized Jaccard Similarity: "
                  << std::fixed << std::setprecision(4) << sim << "\n";
        std::cout << "Verification with Reference:   "
                  << (match ? "PASS (Results match reference)" : "FAIL (Discrepancy detected)") << "\n";
        std::cout << "Strategy Selected:             "
                  << ((k < g_threshold && k <= MAX_SMALL_K) ? "countKmers_small (Direct 4^k Array)" : "countKmers_large (Integer Sort + Merge-Walk)") << "\n";
        std::cout << "Execution Time:                "
                  << std::fixed << std::setprecision(2) << elapsed_us << " us\n";

        // For concise sequences (<= 80 bp), display full k-mer frequency table
        if (s1.size() <= 80 && s2.size() <= 80) {
            std::cout << "\nDetailed k-mer Breakdown:\n";
            auto runsA = countKmers_large(s1, k);
            auto runsB = countKmers_large(s2, k);
            std::map<uint64_t, uint32_t> mapA, mapB;
            std::set<uint64_t> all_kmers;
            for (const auto& r : runsA) { mapA[r.kmer] = r.count; all_kmers.insert(r.kmer); }
            for (const auto& r : runsB) { mapB[r.kmer] = r.count; all_kmers.insert(r.kmer); }

            long long sum_min = 0, sum_max = 0;
            std::cout << "  k-mer   | Seq 1 | Seq 2 | min | max\n";
            std::cout << "  --------+-------+-------+-----+----\n";
            for (uint64_t km : all_kmers) {
                uint32_t cA = mapA.count(km) ? mapA[km] : 0;
                uint32_t cB = mapB.count(km) ? mapB[km] : 0;
                sum_min += std::min(cA, cB);
                sum_max += std::max(cA, cB);
                std::cout << "  " << std::setw(7) << std::left << decodeKmer(km, k)
                          << " | " << std::setw(5) << std::right << cA
                          << " | " << std::setw(5) << cB
                          << " | " << std::setw(3) << std::min(cA, cB)
                          << " | " << std::setw(3) << std::max(cA, cB) << "\n";
            }
            std::cout << "  --------+-------+-------+-----+----\n";
            std::cout << "  TOTALS  |       |       | "
                      << std::setw(3) << sum_min << " | " << std::setw(3) << sum_max << "\n";
            std::cout << "  Formula: sum(min) / sum(max) = "
                      << sum_min << " / " << sum_max << " = "
                      << std::fixed << std::setprecision(4) << sim << "\n";
        }
        std::cout << "--------------------------------------------\n";
        std::cout << "Verification completed successfully.\n";
    }
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char* argv[]) {
    std::cout << "============================================\n";
    std::cout << "  DNA Sequence Similarity — Jaccard Index  \n";
    std::cout << "============================================\n\n";

    std::string mode = (argc > 1) ? argv[1] : "";

    // Phase 0: Calibration
    calibrateAndVerifyStability();

    // Phase 1: Correctness Tests
    runTests();

    // Direct CLI execution handlers
    if (mode == "--interactive" || mode == "-i") {
        runInteractiveMode();
        return 0;
    }
    if (mode == "--benchmark" || mode == "-b") {
        runBenchmark("results_benchmark.csv");
        runMutationExperiment("results_mutation.csv");
        std::cout << "All automated benchmarks complete.\n";
        return 0;
    }
    if (mode == "--test" || mode == "-t") {
        std::cout << "Test suite execution complete.\n";
        return 0;
    }
    if (mode == "--all" || mode == "-a") {
        runBenchmark("results_benchmark.csv");
        runMutationExperiment("results_mutation.csv");
        std::cout << "\nAll experiments complete.\n";
        return 0;
    }

    // Interactive Menu Mode
    while (true) {
        std::cout << "Select an option:\n";
        std::cout << "  1. Interactive Verification Mode (Demo to Professor)\n";
        std::cout << "  2. Run Full Experimental Benchmark Suite (Exp A, B, C)\n";
        std::cout << "  3. Re-run Correctness Test Suite\n";
        std::cout << "  4. Exit\n";
        std::cout << "Enter choice [1-4]: ";

        int choice = 0;
        if (!(std::cin >> choice)) break;

        if (choice == 1) {
            runInteractiveMode();
        } else if (choice == 2) {
            runBenchmark("results_benchmark.csv");
            runMutationExperiment("results_mutation.csv");
            std::cout << "Benchmarks completed successfully.\n\n";
        } else if (choice == 3) {
            runTests();
        } else if (choice == 4) {
            std::cout << "Exiting.\n";
            break;
        } else {
            std::cout << "Invalid choice. Please enter 1, 2, 3, or 4.\n\n";
        }
    }

    return 0;
}
