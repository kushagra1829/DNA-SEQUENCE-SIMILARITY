// ============================================================================
//  DNA Sequence Similarity — Generalized Jaccard Index
//  Performance Programming Assignment
//
//  Design:
//    Two separate k-mer counting strategies are provided:
//      countKmers_small(seq, k)  — optimised for SMALL k (hash-map on
//                                  std::string keys, bucket-friendly for
//                                  short, highly-repeated keys).
//      countKmers_large(seq, k)  — optimised for LARGE k (uses string_view
//                                  for zero-copy key lookup before inserting,
//                                  reduces allocation pressure for long, mostly-
//                                  unique keys that appear infrequently).
//
//    The THRESHOLD between "small" and "large" is determined experimentally
//    at programme start by a calibration phase that measures real execution
//    time on the current machine.  It is NOT a fixed theoretical constant.
//
//  Build (GCC / MinGW, release):
//    g++ -O2 -std=c++17 -o dna_sim main.cpp
//
//  Output files written to the current working directory:
//    results_benchmark.csv    — Experiment A & B data
//    results_mutation.csv     — Experiment C data
//    calibration_log.csv      — raw calibration timings used to pick threshold
// ============================================================================

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------
using FreqMap   = std::unordered_map<std::string, int>;
using Clock     = std::chrono::steady_clock;
using Seconds   = std::chrono::duration<double>;

// ---------------------------------------------------------------------------
// Utility: high-resolution elapsed time
// ---------------------------------------------------------------------------
inline double elapsed_sec(Clock::time_point start) {
    return Seconds(Clock::now() - start).count();
}

// ===========================================================================
// SECTION 1 — k-mer counting  (TWO SEPARATE IMPLEMENTATIONS)
// ===========================================================================

// ---------------------------------------------------------------------------
// Strategy A — SMALL k
//   For small k the number of distinct k-mers is bounded by min(4^k, n-k+1).
//   At k≤8 there are at most 65 536 distinct keys.  Strings are short, hash
//   cost is low, and bucket contention is manageable.  We use a pre-reserved
//   unordered_map and insert std::string keys directly.
// ---------------------------------------------------------------------------
FreqMap countKmers_small(const std::string& seq, int k) {
    FreqMap freq;
    if (k <= 0 || static_cast<int>(seq.size()) < k) return freq;

    // Reserve a reasonable number of buckets to avoid rehashing.
    // For small k (k <= 12), 4^k is the maximum number of distinct k-mers.
    // For k >= 13, 4^k exceeds practical DNA distinct k-mers so we bound by n-k+1.
    size_t n_kmers = static_cast<size_t>(seq.size() - k + 1);
    size_t reserve_count = n_kmers;
    if (k <= 12) {
        size_t p4 = 1ULL << (2 * k); // 4^k = 2^(2k)
        if (p4 < reserve_count) reserve_count = p4;
    }
    freq.reserve(reserve_count);

    const int n = static_cast<int>(seq.size());
    for (int i = 0; i <= n - k; ++i) {
        freq[seq.substr(i, k)]++;
    }
    return freq;
}

// ---------------------------------------------------------------------------
// Strategy B — LARGE k
//   For large k, keys are long and repeated k-mers waste heap memory if
//   materialized as std::string prior to map insertion.  Here we use zero-copy
//   std::string_view stored in a contiguous std::vector, sort lexicographically
//   with high cache locality, and perform run-length consolidation.
//   std::string keys are instantiated ONLY ONCE per UNIQUE k-mer, eliminating
//   heap churn and fragmentation for repeated sequences.
// ---------------------------------------------------------------------------
FreqMap countKmers_large(const std::string& seq, int k) {
    FreqMap freq;
    if (k <= 0 || static_cast<int>(seq.size()) < k) return freq;

    const int n = static_cast<int>(seq.size());
    const int num_kmers = n - k + 1;

    std::vector<std::string_view> views;
    views.reserve(num_kmers);
    for (int i = 0; i < num_kmers; ++i) {
        views.emplace_back(seq.data() + i, k);
    }

    std::sort(views.begin(), views.end());

    // Single pass run-length consolidation: allocate std::string only for unique k-mers
    freq.reserve(num_kmers);
    int run_len = 1;
    for (size_t i = 1; i < views.size(); ++i) {
        if (views[i] == views[i - 1]) {
            ++run_len;
        } else {
            freq.emplace(std::string(views[i - 1]), run_len);
            run_len = 1;
        }
    }
    if (!views.empty()) {
        freq.emplace(std::string(views.back()), run_len);
    }
    return freq;
}

// ===========================================================================
// SECTION 2 — Calibration: find the experimental threshold
// ===========================================================================

// We benchmark countKmers_small vs countKmers_large on a shared fixed test
// sequence for k = 1..32, record average timings, and return the first k
// at which countKmers_large becomes faster (or stays within 5% of small).
// If large is never faster, we default to the midpoint k=16 as a safe fallback.

struct CalibrationResult {
    int threshold;                       // k value at which we switch
    std::vector<double> time_small_ms;   // avg time (ms) per k for small
    std::vector<double> time_large_ms;   // avg time (ms) per k for large
};

CalibrationResult calibrate(int calib_len = 50000, int repeats = 5) {
    std::cout << "[Calibration] Generating test sequence of length "
              << calib_len << " (repeats=" << repeats << ")...\n";

    // Fixed seed for reproducibility
    std::mt19937 rng(42);
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<int> base_dist(0, 3);

    std::string seq(calib_len, 'A');
    for (char& c : seq) c = bases[base_dist(rng)];

    CalibrationResult cr;
    cr.time_small_ms.resize(33, 0.0);
    cr.time_large_ms.resize(33, 0.0);

    std::ofstream log("calibration_log.csv");
    log << "k,strategy,avg_time_ms\n";

    for (int k = 1; k <= 32; ++k) {
        // --- small strategy ---
        double total_s = 0.0;
        bool small_ok = true;
        for (int r = 0; r < repeats; ++r) {
            try {
                auto t0 = Clock::now();
                auto fm = countKmers_small(seq, k);
                total_s += elapsed_sec(t0);
                (void)fm.size();
            } catch (const std::bad_alloc&) {
                small_ok = false; break;
            }
        }
        cr.time_small_ms[k] = small_ok ? (total_s / repeats) * 1000.0 : -1.0;

        // --- large strategy ---
        double total_l = 0.0;
        bool large_ok = true;
        for (int r = 0; r < repeats; ++r) {
            try {
                auto t0 = Clock::now();
                auto fm = countKmers_large(seq, k);
                total_l += elapsed_sec(t0);
                (void)fm.size();
            } catch (const std::bad_alloc&) {
                large_ok = false; break;
            }
        }
        cr.time_large_ms[k] = large_ok ? (total_l / repeats) * 1000.0 : -1.0;

        log << k << ",small," << cr.time_small_ms[k] << "\n";
        log << k << ",large," << cr.time_large_ms[k] << "\n";

        std::string s_str = small_ok
            ? (std::to_string(static_cast<int>(cr.time_small_ms[k]*1000)/1000.0).substr(0,6) + "ms")
            : "OOM";
        std::string l_str = large_ok
            ? (std::to_string(static_cast<int>(cr.time_large_ms[k]*1000)/1000.0).substr(0,6) + "ms")
            : "OOM";

        std::cout << "  k=" << std::setw(2) << k
                  << "  small=" << s_str
                  << "  large=" << l_str << "\n";
    }

    // Find first k where large is faster (within 5% margin to avoid noise)
    // If no crossover is found, keep using small for all k (threshold = 33).
    cr.threshold = 33; // default: small strategy wins across all k=1..32
    for (int k = 1; k <= 32; ++k) {
        double small_t = cr.time_small_ms[k];
        double large_t = cr.time_large_ms[k];
        if (small_t > 0 && large_t > 0 && large_t < small_t * 0.95) {   // large is at least 5% faster
            cr.threshold = k;
            break;
        }
    }

    if (cr.threshold > 32) {
        std::cout << "[Calibration] No crossover found: small strategy is faster"
                     " for ALL k=1..32 on this machine.\n";
        std::cout << "[Calibration] Threshold set to 33 (small used for all k).\n\n";
    } else {
        std::cout << "[Calibration] Experimental threshold = k=" << cr.threshold
                  << "  (large strategy used for k>=" << cr.threshold << ")\n\n";
    }
    return cr;
}

// ===========================================================================
// SECTION 3 — Unified dispatcher
// ===========================================================================

// Global threshold set after calibration
int g_threshold = 17;

FreqMap countKmers(const std::string& seq, int k) {
    if (k < g_threshold) {
        return countKmers_small(seq, k);
    } else {
        return countKmers_large(seq, k);
    }
}

// ===========================================================================
// SECTION 4 — Generalized Jaccard similarity
// ===========================================================================

double generalizedJaccard(const std::string& seqA,
                           const std::string& seqB,
                           int k) {
    if (k <= 0) return 0.0;

    FreqMap fA = countKmers(seqA, k);
    FreqMap fB = countKmers(seqB, k);

    if (fA.empty() && fB.empty()) return 1.0; // both too short: identical edge case
    if (fA.empty() || fB.empty()) return 0.0;

    long long sum_min = 0, sum_max = 0;

    // Iterate over all k-mers in A
    for (auto& [kmer, cA] : fA) {
        auto it = fB.find(kmer);
        int cB = (it != fB.end()) ? it->second : 0;
        sum_min += std::min(cA, cB);
        sum_max += std::max(cA, cB);
    }
    // Add k-mers present only in B
    for (auto& [kmer, cB] : fB) {
        if (fA.find(kmer) == fA.end()) {
            sum_max += cB; // min is 0
        }
    }

    if (sum_max == 0) return 1.0; // degenerate: both empty after filtering
    return static_cast<double>(sum_min) / static_cast<double>(sum_max);
}

// Timed variant: returns {similarity, time_count_ms, time_jaccard_ms, total_ms}
struct TimedResult {
    double similarity;
    double time_count_ms;   // time to build both frequency maps
    double time_jaccard_ms; // time to compute Jaccard from existing maps
    double total_ms;
};

TimedResult timedJaccard(const std::string& seqA,
                          const std::string& seqB,
                          int k) {
    if (k <= 0) return {0.0, 0.0, 0.0, 0.0};

    // --- Phase 1: count k-mers ---
    auto t0 = Clock::now();
    FreqMap fA = countKmers(seqA, k);
    FreqMap fB = countKmers(seqB, k);
    double count_ms = elapsed_sec(t0) * 1000.0;

    // --- Phase 2: compute Jaccard ---
    auto t1 = Clock::now();
    long long sum_min = 0, sum_max = 0;

    if (fA.empty() && fB.empty()) {
        return {1.0, count_ms, 0.0, count_ms};
    }
    for (auto& [kmer, cA] : fA) {
        auto it = fB.find(kmer);
        int cB = (it != fB.end()) ? it->second : 0;
        sum_min += std::min(cA, cB);
        sum_max += std::max(cA, cB);
    }
    for (auto& [kmer, cB] : fB) {
        if (fA.find(kmer) == fA.end()) {
            sum_max += cB;
        }
    }

    double sim = (sum_max == 0) ? 1.0
                                : static_cast<double>(sum_min) / sum_max;
    double jaccard_ms = elapsed_sec(t1) * 1000.0;

    return {sim, count_ms, jaccard_ms, count_ms + jaccard_ms};
}

// ===========================================================================
// SECTION 5 — DNA sequence generation and mutation
// ===========================================================================

std::string generateDNA(int length, std::mt19937& rng) {
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<int> dist(0, 3);
    std::string seq(length, 'A');
    for (char& c : seq) c = bases[dist(rng)];
    return seq;
}

// Introduce exactly `num_mutations` point mutations (substitutions only).
// Returns a new string; does not modify the original.
std::string mutate(const std::string& seq, int num_mutations, std::mt19937& rng) {
    std::string mutated = seq;
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<int> pos_dist(0, static_cast<int>(seq.size()) - 1);
    std::uniform_int_distribution<int> base_dist(0, 3);

    int applied = 0;
    while (applied < num_mutations) {
        int pos  = pos_dist(rng);
        char nb  = bases[base_dist(rng)];
        if (nb != mutated[pos]) {         // must be a true substitution
            mutated[pos] = nb;
            ++applied;
        }
    }
    return mutated;
}

// ===========================================================================
// SECTION 6 — Correctness test cases
// ===========================================================================

void runTests() {
    std::cout << "=== Correctness Tests ===\n\n";

    // Test 1: identical sequences → similarity must be 1.0
    {
        std::string s = "ATTCGTA";
        double sim = generalizedJaccard(s, s, 3);
        std::cout << "[T1] Identical sequences, k=3: sim=" << sim
                  << (std::abs(sim - 1.0) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 2: completely different alphabet-disjoint sequences
    //   All-A vs All-C: k-mer sets are disjoint → 0
    {
        std::string a(10, 'A'), c(10, 'C');
        double sim = generalizedJaccard(a, c, 3);
        std::cout << "[T2] Completely different (AAAA vs CCCC), k=3: sim=" << sim
                  << (sim == 0.0 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 3: hand-calculated example
    //   seqA = "AACT", k=2 → AA:1, AC:1, CT:1
    //   seqB = "AACT", k=2 → same → Jaccard = 1.0
    {
        std::string a = "AACT", b = "AACT";
        double sim = generalizedJaccard(a, b, 2);
        std::cout << "[T3] AACT vs AACT, k=2: sim=" << sim
                  << (std::abs(sim - 1.0) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 4: hand-calculated example with partial overlap
    //   seqA = "AACG", k=2 → {AA:1, AC:1, CG:1}
    //   seqB = "AATG", k=2 → {AA:1, AT:1, TG:1}
    //   Intersection: AA → min(1,1)=1   sum_min=1
    //   Union:        AA→max=1, AC→1, CG→1, AT→1, TG→1  sum_max=5
    //   Expected Jaccard = 1/5 = 0.2
    {
        std::string a = "AACG", b = "AATG";
        double sim = generalizedJaccard(a, b, 2);
        std::cout << "[T4] AACG vs AATG, k=2: sim=" << sim
                  << (std::abs(sim - 0.2) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 5: repeated pattern
    //   seqA = "ATATAT", k=2 → {AT:3, TA:2}  (AT appears at 0,2,4; TA at 1,3)
    //   seqB = "ATATAT" same → Jaccard = 1.0
    {
        std::string s = "ATATAT";
        double sim = generalizedJaccard(s, s, 2);
        std::cout << "[T5] Repeated pattern ATATAT vs itself, k=2: sim=" << sim
                  << (std::abs(sim - 1.0) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 6: k larger than sequence → both maps empty → sim=1 by convention
    {
        std::string s = "ACG";
        double sim = generalizedJaccard(s, s, 10);
        std::cout << "[T6] k > seq length (ACG, k=10): sim=" << sim
                  << (std::abs(sim - 1.0) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 7: k=1, identical distribution
    //   "ACGT" vs "TGCA": each base appears exactly once → Jaccard = 1
    {
        std::string a = "ACGT", b = "TGCA";
        double sim = generalizedJaccard(a, b, 1);
        std::cout << "[T7] ACGT vs TGCA, k=1 (same base freq): sim=" << sim
                  << (std::abs(sim - 1.0) < 1e-9 ? "  PASS" : "  FAIL") << "\n";
    }

    // Test 8: k=32 on medium sequence (stress)
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
// SECTION 7 — Experiment A & B: benchmark sequence length vs k
// ===========================================================================

void runBenchmark(const std::string& csv_path) {
    // Sequence lengths to test.
    // 1,000,000 bp requires ~2 GB RAM for k-mer maps at large k; it is
    // attempted but skipped gracefully if the system runs out of memory.
    const std::vector<int> lengths = {1000, 10000, 100000, 500000, 1000000};
    // k values to test
    std::vector<int> k_values;
    for (int k = 1; k <= 32; ++k) k_values.push_back(k);

    const int REPEATS = 5;  // number of timing repetitions per cell

    std::ofstream csv(csv_path);
    csv << "seq_len,k,strategy,repetition,"
        << "similarity,time_count_ms,time_jaccard_ms,total_ms\n";

    std::mt19937 rng(1234);

    std::cout << "=== Experiment A & B: Benchmark ===\n";

    for (int len : lengths) {
        std::cout << "\n  Generating sequences of length " << len << "...\n";

        std::string seqA, seqB;
        try {
            seqA = generateDNA(len, rng);
            seqB = generateDNA(len, rng);
        } catch (const std::bad_alloc&) {
            std::cout << "  [SKIP] Out of memory generating sequences of length "
                      << len << ". Skipping.\n";
            continue;
        }

        bool oom_hit = false;
        for (int k : k_values) {
            if (oom_hit) break;
            std::string strategy = (k < g_threshold) ? "small" : "large";

            // Warm-up run (not recorded)
            try { auto r = timedJaccard(seqA, seqB, k); (void)r; }
            catch (const std::bad_alloc&) { oom_hit = true; break; }

            for (int rep = 1; rep <= REPEATS; ++rep) {
                try {
                    auto r = timedJaccard(seqA, seqB, k);
                    csv << len << "," << k << "," << strategy << "," << rep << ","
                        << std::fixed << std::setprecision(6)
                        << r.similarity << ","
                        << r.time_count_ms << ","
                        << r.time_jaccard_ms << ","
                        << r.total_ms << "\n";
                } catch (const std::bad_alloc&) {
                    std::cout << "  [OOM] Out of memory at len=" << len
                              << " k=" << k << ". Skipping remaining k for this length.\n";
                    oom_hit = true;
                    break;
                }
            }
            if (!oom_hit)
                std::cout << "    k=" << std::setw(2) << k
                          << " [" << strategy << "]  done\n";
        }
    }

    std::cout << "\n[Benchmark] Results written to: " << csv_path << "\n\n";
}

// ===========================================================================
// SECTION 8 — Experiment C: mutation sensitivity
// ===========================================================================

void runMutationExperiment(const std::string& csv_path) {
    const int SEQ_LEN = 100000;
    const std::vector<int> k_values = {3, 8, 16, 24, 32};
    // Mutation percentages (0%, 0.1%, 0.5%, 1%, 2%, 5%, 10%, 20%, 50%)
    const std::vector<double> mut_pcts = {
        0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0
    };
    const int REPEATS = 3;

    std::cout << "=== Experiment C: Mutation Sensitivity ===\n";
    std::cout << "  Sequence length: " << SEQ_LEN << "\n\n";

    std::mt19937 rng(9999);
    std::string original = generateDNA(SEQ_LEN, rng);

    std::ofstream csv(csv_path);
    csv << "seq_len,k,mutation_pct,num_mutations,repetition,"
        << "similarity,total_ms\n";

    for (int k : k_values) {
        std::cout << "  k=" << k << "\n";
        for (double pct : mut_pcts) {
            int num_mut = static_cast<int>(std::round(SEQ_LEN * pct / 100.0));
            num_mut = std::min(num_mut, SEQ_LEN); // safety cap

            for (int rep = 1; rep <= REPEATS; ++rep) {
                std::string mutated = (num_mut == 0)
                    ? original
                    : mutate(original, num_mut, rng);

                auto r = timedJaccard(original, mutated, k);

                csv << SEQ_LEN << "," << k << ","
                    << pct << "," << num_mut << "," << rep << ","
                    << std::fixed << std::setprecision(6)
                    << r.similarity << "," << r.total_ms << "\n";
            }
            std::cout << "    mut=" << std::setw(5) << pct << "%"
                      << "  num_mut=" << std::setw(6) << num_mut << "\n";
        }
    }

    std::cout << "\n[Mutation] Results written to: " << csv_path << "\n\n";
}

// ===========================================================================
// SECTION 9 — Interactive Verification Mode & Reference Implementation
// ===========================================================================

// Straightforward reference implementation using standard sorted map
// for cross-verification of optimized implementations.
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

        // Compute similarity with timing
        auto t0 = Clock::now();
        double sim = generalizedJaccard(s1, s2, k);
        double elapsed_us = elapsed_sec(t0) * 1e6;

        // Cross-verify with reference implementation
        double ref_sim = referenceJaccard(s1, s2, k);
        bool match = (std::abs(sim - ref_sim) < 1e-9);

        std::cout << "\n--------------------------------------------\n";
        std::cout << "Generalized Jaccard Similarity: "
                  << std::fixed << std::setprecision(4) << sim << "\n";
        std::cout << "Verification with Reference:   "
                  << (match ? "PASS (Results match reference)" : "FAIL (Discrepancy detected)") << "\n";
        std::cout << "Strategy Selected:             "
                  << ((k < g_threshold) ? "countKmers_small (Hash Table)" : "countKmers_large (string_view Sort)") << "\n";
        std::cout << "Execution Time:                "
                  << std::fixed << std::setprecision(2) << elapsed_us << " us\n";

        // For concise sequences (<= 80 bp), display full k-mer frequency table
        if (s1.size() <= 80 && s2.size() <= 80) {
            std::cout << "\nDetailed k-mer Breakdown:\n";
            FreqMap fA = countKmers(s1, k);
            FreqMap fB = countKmers(s2, k);
            std::set<std::string> all_kmers;
            for (const auto& [km, _] : fA) all_kmers.insert(km);
            for (const auto& [km, _] : fB) all_kmers.insert(km);

            long long sum_min = 0, sum_max = 0;
            std::cout << "  k-mer   | Seq 1 | Seq 2 | min | max\n";
            std::cout << "  --------+-------+-------+-----+----\n";
            for (const auto& km : all_kmers) {
                int cA = fA.count(km) ? fA[km] : 0;
                int cB = fB.count(km) ? fB[km] : 0;
                sum_min += std::min(cA, cB);
                sum_max += std::max(cA, cB);
                std::cout << "  " << std::setw(7) << std::left << km
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

    // Check CLI argument shortcuts
    std::string mode = "";
    if (argc > 1) {
        mode = argv[1];
    }

    // ------------------------------------------------------------------
    // Phase 0 — Experimental calibration to determine strategy threshold
    // ------------------------------------------------------------------
    CalibrationResult cal = calibrate(30000, 3);
    g_threshold = cal.threshold;

    std::cout << ">>> Using threshold k=" << g_threshold
              << " (small strategy: k<" << g_threshold
              << ", large strategy: k>=" << g_threshold << ")\n\n";

    // ------------------------------------------------------------------
    // Phase 1 — Correctness Verification Tests
    // ------------------------------------------------------------------
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
        if (!(std::cin >> choice)) {
            break;
        }

        if (choice == 1) {
            runInteractiveMode();
        } else if (choice == 2) {
            runBenchmark("results_benchmark.csv");
            runMutationExperiment("results_mutation.csv");
            std::cout << "Benchmarks completed successfully.\n\n";
        } else if (choice == 3) {
            runTests();
        } else if (choice == 4) {
            std::cout << "Exiting DNA Sequence Similarity program.\n";
            break;
        } else {
            std::cout << "Invalid choice. Please enter 1, 2, 3, or 4.\n\n";
        }
    }

    return 0;
}

