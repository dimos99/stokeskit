#include "YC.h"
#include "cache_hash.h"

#include "common.h"
#include <fstream>
#include <filesystem>
#include <sstream>
#include <iomanip>

// Add these helper functions before Cache_YC struct
namespace {
    // Maximum values for n, p, q to determine vector sizes
    constexpr int MAX_N = 600;
    constexpr int MAX_P = 600;
    constexpr int MAX_Q = 600;
    
    // Convert (n,p,q) to single index
    inline size_t npq_to_index(int n, int p, int q) {
        return ((size_t)n * MAX_P * MAX_Q) + ((size_t)p * MAX_Q) + (size_t)q;
    }

    // Get vector size needed for cache
    constexpr size_t get_cache_size() {
        return (size_t)MAX_N * MAX_P * MAX_Q;
    }

    // Cache file paths
    const std::string CACHE_DIR = ".cache";
    const std::string P_CACHE_FILE = CACHE_DIR + "/P_YC_cache.bin";
    const std::string V_CACHE_FILE = CACHE_DIR + "/V_YC_cache.bin";
    const std::string Q_CACHE_FILE = CACHE_DIR + "/Q_YC_cache.bin";
}

/**
 * @struct GCoeffs_YC
 * @brief Calculates g1, g2 and g3 coefficients for the YC scalar resistance functions.
 *
 * These calculations follow the formulas in the paper by Jeffrey and Onishi (1984).
 * They are intermediate calculations used in the YC11 and YC12 functions.
 *
 * @param l Radius ratio parameter used in calculations
 *
 * Members:
 * @var g2 Second coefficient
 * @var g3 Third coefficient
 */
struct GCoeffs_YC {
    double g2, g3, g4, g5;
    GCoeffs_YC(double l) {
        // Precompute common terms
        const double l2 = l * l;
        const double one_plus_l = 1.0 + l;
        const double l_plus_1_inv = 1.0 / one_plus_l;
        const double l_plus_1_inv_to_the_4 = l_plus_1_inv * l_plus_1_inv * l_plus_1_inv * l_plus_1_inv;
        
        // Calculate coefficients using precomputed terms
        g2 = 0.4 * l * l_plus_1_inv;
        g3 = 0.008 * (8.0 + 6.0 * l + 33.0 * l2) * l_plus_1_inv;
        g4 = 0.8 * l2 * l_plus_1_inv_to_the_4;
        g5 = 0.016 * l * (43.0 - 24.0 * l + 43.0 * l2) * l_plus_1_inv_to_the_4; // TODO: Check this
    }
};

// Global cache structure for all calculations
/**
 * @brief Singleton class for global caching of frequently used calculations
 * 
 * This structure implements a singleton pattern to provide global caching functionality
 * for various mathematical computations. It helps avoid redundant calculations by
 * storing previously computed values.
 * 
 * @details The cache contains the following components:
 * - P_cache: Stores P-function values indexed by three integers
 * - V_cache: Stores V-function values indexed by three integers
 * - fk_cache: Stores f_k function values indexed by integer and double pair
 * - CY_cache: Combined cache for C11Y and C12Y functions, indexed by type flag and l parameter
 * - g_coeffs_cache: Stores G-coefficients indexed by double
 * - comb_cache: Stores combinatorial values indexed by integer pairs
 * 
 * The singleton pattern ensures that only one instance of the cache exists throughout
 * the program's lifetime, providing a centralized storage for computed values.
 * 
 * Usage:
 * @code
 * Cache_YC& cache = Cache_YC::getInstance();
 * cache.clear(); // To clear all cached values
 * @endcode
 */
struct Cache_YC {
    static const size_t MAX_CACHE_SIZE = 100000000000;

    std::vector<std::pair<bool, double>> P_cache; // (valid, value) pairs
    std::vector<std::pair<bool, double>> V_cache;
    std::vector<std::pair<bool, double>> Q_cache;
    std::unordered_map<std::pair<int, double>, double, pair_int_double_hash> fk_cache; // Keep this as unordered_map
    std::unordered_map<std::pair<bool, double>, double, pair_bool_double_hash> CY_cache; // Keep this as unordered_map
    std::unordered_map<double, GCoeffs_YC> g_coeffs_cache; // Keep this as unordered_map
    
    Cache_YC() : 
        P_cache(get_cache_size(), {false, 0.0}),
        V_cache(get_cache_size(), {false, 0.0}),
        Q_cache(get_cache_size(), {false, 0.0}),
        fk_cache(),
        CY_cache(),
        g_coeffs_cache()
    {
        cache_utils::ensure_cache_dir(CACHE_DIR);
        auto& stats = cache_utils::get_stats();
        stats.cache_modified = false;
        
        cache_utils::read_cache_file(P_CACHE_FILE, P_cache, MAX_N, MAX_P, MAX_Q);
        cache_utils::read_cache_file(V_CACHE_FILE, V_cache, MAX_N, MAX_P, MAX_Q);
        cache_utils::read_cache_file(Q_CACHE_FILE, Q_cache, MAX_N, MAX_P, MAX_Q);
#ifdef ENABLE_PROFILING
        updateMemoryStats();
#endif
    }

    static Cache_YC& getInstance() {
        static Cache_YC instance;
        return instance;
    }

    void clear() {
        std::fill(P_cache.begin(), P_cache.end(), std::make_pair(false, 0.0));
        std::fill(V_cache.begin(), V_cache.end(), std::make_pair(false, 0.0));
        std::fill(Q_cache.begin(), Q_cache.end(), std::make_pair(false, 0.0));
        fk_cache.clear();
        g_coeffs_cache.clear();
#ifdef ENABLE_PROFILING
        updateMemoryStats();
#endif
    }

    // Add template method for inserting into unordered_map caches
    template<typename K, typename V, typename H>
    void insertWithCheck(std::unordered_map<K,V,H>& cache, const K& key, const V& value) {
        cache[key] = value;
        if (cache.size() > MAX_CACHE_SIZE) {
            std::cerr << "Warning: YC cache size exceeded maximum limit. Clearing cache..." << std::endl;
            auto it = cache.begin();
            size_t to_remove = cache.size() / 2;
            for (size_t i = 0; i < to_remove && it != cache.end(); ++i) {
                it = cache.erase(it);
            }
        }
#ifdef ENABLE_PROFILING
        updateMemoryStats();
#endif
    }

#ifdef ENABLE_PROFILING
    void updateMemoryStats() {
        size_t total_bytes = 0;
        size_t total_entries = 0;

        // P_cache memory
        total_bytes += sizeof(bool) * P_cache.size();
        total_bytes += sizeof(double) * P_cache.size();
        total_entries += std::count_if(P_cache.begin(), P_cache.end(), 
                                     [](const auto& p) { return p.first; });

        // V_cache memory
        total_bytes += sizeof(bool) * V_cache.size();
        total_bytes += sizeof(double) * V_cache.size();
        total_entries += std::count_if(V_cache.begin(), V_cache.end(), 
                                     [](const auto& p) { return p.first; });

        // Q_cache memory
        total_bytes += sizeof(bool) * Q_cache.size();
        total_bytes += sizeof(double) * Q_cache.size();
        total_entries += std::count_if(Q_cache.begin(), Q_cache.end(), 
                                     [](const auto& p) { return p.first; });

        MemoryProfiler::record("Cache_YC", total_bytes, total_entries);
    }
#endif

    ~Cache_YC() {
        auto& stats = cache_utils::get_stats();
        if (stats.cache_modified) {
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    cache_utils::write_cache_file(P_CACHE_FILE, P_cache, MAX_N, MAX_P, MAX_Q);
                }
                #pragma omp section
                {
                    cache_utils::write_cache_file(V_CACHE_FILE, V_cache, MAX_N, MAX_P, MAX_Q);
                }
                #pragma omp section
                {
                    cache_utils::write_cache_file(Q_CACHE_FILE, Q_cache, MAX_N, MAX_P, MAX_Q);
                }
            }
            cache_utils::print_io_stats();
        }
    }

private:
    Cache_YC(const Cache_YC&) = delete;
    Cache_YC& operator=(const Cache_YC&) = delete;
};


/**
 * @brief Calculates the pure P_{n,p,q} function for the YC scalar function in a two-body resistance problem.
 * 
 * 
 * @details
 * This function implements a recursive formula for calculating P_{n,p,q} coefficients
 * used in the series expansion of the two-body problem solutions. It uses memoization
 * through a global cache to optimize recursive calls. See equation (3.6, 3.9) in Jeffrey and Onishi (1984).
 * 
 * 
 * @param n Integer parameter n in P_{n,p,q}
 * @param p Integer parameter p in P_{n,p,q}
 * @param q Integer parameter q in P_{n,p,q}
 * 
 * @return The value of P_{n,p,q} for the given parameters
 * 
 * @note The function uses caching to avoid redundant calculations
 * @note Base case is when p = q = 0, where P_{n,0,0} = 1 if n = 1, and 0 otherwise
 * 
 * @see V_npq_YC Function used in the recursive formula
 * @see Q_npq_YC Function used in the recursive formula
 * @see Cache_YC::P_cache Cache storage for memoization
 */
double P_npq_YC(int n, int p, int q) {

        if (n < 0 || p < 0 || q < 0) {
        return 0.0;
    }


#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    if (n >= MAX_N || p >= MAX_P || q >= MAX_Q) {
        throw std::runtime_error("Cache indices out of bounds in P_npq_YC");
    }

    size_t idx = npq_to_index(n, p, q);
    auto& cache_entry = Cache_YC::getInstance().P_cache[idx];
    
    if(cache_entry.first) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("P_npq_YC", duration, cache_hit);
#endif
        return cache_entry.second;
    }

    // Base case
    if(p == 0 && q == 0) {
        double P = 0.0;
        cache_entry.first = true;
        cache_entry.second = P;
        cache_utils::get_stats().cache_modified = true;  // Add this line
#ifdef ENABLE_PROFILING
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("P_npq_YC", duration, cache_hit);
#endif
        return P;
    }

    double P = 0.0;
    double pre_term1 = (2.0*n + 1.0) / (2.0*(n + 1.0));
    double pre_term2 = (n* (2.0*n - 1.0)) / (2.0*(n + 1.0));
    double pre_term3 = (n*(4.0*n*n - 1.0)) / (2.0*(n + 1.0));
    double pre_term4 = - (2.0*(4.0*n * n - 1.0 )) / (3.0*(n + 1.0));
    for(int s = 1; s <= q; ++s) {
        double term1 = pre_term1 * (3.0 * (n + s) - (n * s + 1.0) * (2.0 * n * s - s - n + 2.0)) /
                                    (s * (n + s) * (2.0 * s - 1.0));
        double term2 = pre_term2;
        double term3 = pre_term3 * 1.0 / (2.0 * s + 1.0);
        double term4 = pre_term4;
        P += comb(n + s, n + 1) * (
            term1 * P_npq_YC(s, q - s, p - n + 1) +
            term2 * P_npq_YC(s, q - s, p - n - 1) +
            term3 * V_npq_YC(s, q - s - 2, p - n + 1) +
            term4 * Q_npq_YC(s, q - s - 1, p - n + 1)
        );
    }

    cache_entry.first = true;
    cache_entry.second = P;
    cache_utils::get_stats().cache_modified = true;  // Add this line
#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("P_npq_YC", duration, cache_hit);
#endif
    return P;
}



/**
 * @brief Calculates the V_{n,p,q} pure function for the YC scalar function in a two-body resistance problem.
 * 
 * @details
 * This function implements a recursive formula for computing V_{n,p,q} coefficients,
 * which are used in the series expansion of the two-body problem.
 * It uses memoization through a global cache to improve performance.
 * 
 * @param n Integer parameter n in V_{n,p,q}
 * @param p Integer parameter p in V_{n,p,q}
 * @param q Integer parameter q in V_{n,p,q}
 * 
 * @return The value of the V_{n,p,q} coefficient
 * 
 * @note The function uses Cache_YC for memoization
 * @note Depends on P_npq_YC function and comb (combination) function
 * @note This function is wrong in Jeffrey & Onishi (1984) paper, the correct formula is in Townsend (2023)
 * 
 * @see P_npq_YC(), comb()
 */
double V_npq_YC(int n, int p, int q) {

    if (n < 0 || p < 0 || q < 0) {
        return 0.0;
    }



#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    if (n >= MAX_N || p >= MAX_P || q >= MAX_Q) {
        throw std::runtime_error("Cache indices out of bounds in V_npq_YC");
    }

    size_t idx = npq_to_index(n, p, q);
    auto& cache_entry = Cache_YC::getInstance().V_cache[idx];
    
    if(cache_entry.first) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("V_npq_YC", duration, cache_hit);
#endif
        return cache_entry.second;
    }

    if(p == 0 && q == 0) {
        double V = 0.0;
        cache_entry.first = true;
        cache_entry.second = V;
        cache_utils::get_stats().cache_modified = true;  // Add this line
#ifdef ENABLE_PROFILING
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("V_npq_YC", duration, cache_hit);
#endif
        return V;
    }

    double term = 2.0 * n / ((n + 1.0)*(2.0*n + 3.0));
    double V = 0.0;
    for(int s = 1; s <= q; ++s) {
        V += comb(n + s, n + 1) * P_npq_YC(s, q - s, p - n - 1);
    }

    V = P_npq_YC(n, p, q) + term * V;
    cache_entry.first = true;
    cache_entry.second = V;
    cache_utils::get_stats().cache_modified = true;  // Add this line
#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("V_npq_YC", duration, cache_hit);
#endif
    return V;
}


/** 
 * @brief Calculates the Q_{n,p,q} pure function for the YC scalar function in a two-body resistance problem.
 * 
 * @details
 * This function implements a recursive formula for computing Q_{n,p,q} coefficients,
 * which are used in the series expansion of the two-body problem.
 * It uses memoization through a global cache to improve performance.
 * 
 * @param n Integer parameter n in Q_{n,p,q}
 * @param p Integer parameter p in Q_{n,p,q}
 * @param q Integer parameter q in Q_{n,p,q}
 * 
 * @return The value of the Q_{n,p,q} coefficient
 * 
 * @note The function uses Cache_YC for memoization
 * @note Depends on P_npq_YC function and comb (combination) function
 * 
 * @see P_npq_YC(), comb()
 */
double Q_npq_YC(int n, int p, int q) {

    if (n < 0 || p < 0 || q < 0) {
        return 0.0;
    }

#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    if (n >= MAX_N || p >= MAX_P || q >= MAX_Q) {
        throw std::runtime_error("Cache indices out of bounds in Q_npq_YC");
    }

    size_t idx = npq_to_index(n, p, q);
    auto& cache_entry = Cache_YC::getInstance().Q_cache[idx];
    
    if(cache_entry.first) {
        // std::cout << "Cache hit for Q_npq_YC(" << n << ", " << p << ", " << q << ")" << std::endl;
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("Q_npq_YC", duration, cache_hit);
#endif
        return cache_entry.second;
    }

    if(p == 0 && q == 0) {
        double Q = n == 1 ? 1.0 : 0.0;
        cache_entry.first = true;
        cache_entry.second = Q;
        cache_utils::get_stats().cache_modified = true;  // Add this line
#ifdef ENABLE_PROFILING
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("Q_npq_YC", duration, cache_hit);
#endif
        return Q;
    }

    double pre_term1 = 1.0 / (n + 1.0);
    double pre_term2 = - 3.0 / (2.0 * n * (n + 1.0));
    double Q = 0.0;
    for(int s = 1; s <= q; ++s) {
        Q += comb(n + s, n + 1) * (
            pre_term1 * s * Q_npq_YC(s, q - s - 1, p - n) +
            pre_term2 / s * P_npq_YC(s, q - s, p - n)
        );
    }

    cache_entry.first = true;
    cache_entry.second = Q;
    cache_utils::get_stats().cache_modified = true;  // Add this line

#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("Q_npq_YC", duration, cache_hit);
#endif
    return Q;
}


/**
 * @brief Calculates the fk function for the YC scalar function in a two-body resistance problem.
 * 
 * @details
 * This function implements a recursive calculation of fk with caching mechanism
 * for improved performance. It uses equation 3.15 in Jeffrey and Onishi (1984).
 * 
 * @param k Integer parameter representing the order of the function
 * @param l Double parameter representing the ratio of sphere radii
 * 
 * @return Double value of the calculated fk function
 * 
 * @note The function uses a global cache to store previously calculated values
 * 
 * @see P_npq_YC, Cache_YC
 */
double fk_YC(int k, double l) {
#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    auto& cache = Cache_YC::getInstance().fk_cache;
    auto key = std::make_pair(k, l);
    auto it = cache.find(key);
    if(it != cache.end()) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("fk_YC", duration, cache_hit);
#endif
        return it->second;
    }

    double result = 0.0;    


    for(int q = 0; q <= k; ++q) {
        // #ifdef _OPENMP
        // std::cout << "[fk_YC debug] Thread " << omp_get_thread_num() 
        //           << " handling q = " << q 
        //           << std::endl;
        // #endif
        result += Q_npq_YC(1, k - q, q) * std::pow(l, q + k%2);
    }
    result *= std::pow(2.0, k);

    Cache_YC::getInstance().insertWithCheck(cache, key, result);
#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("fk_YC", duration, cache_hit);
#endif
    return result;
}

// ---------------------- YC11 functions ----------------------

// Lubrication expression

/**
 * @brief Computes the C11Y coefficient for the two-body problem.
 * 
 * @details
 * This function calculates the C11Y coefficient using a series expansion with cached values
 * for optimization. It utilizes g-coefficients and implements a convergence check based on
 * relative and absolute tolerances. See equation 18a of Townsend (2023) for the formula.
 *
 * @param l The ratio of the two sphere radii
 * @param maxIter Maximum number of iterations for the series expansion
 * @param rtol Relative tolerance for convergence check
 * @param atol Absolute tolerance for convergence check
 * 
 * @return The computed C11Y coefficient value
 * 
 * @note The function uses caching mechanisms through Cache_YC for both C11Y values
 *       and g-coefficients to improve performance on repeated calls
 * @warning May print a warning to stderr if convergence is not reached within maxIter iterations
 * 
 * @see GCoeffs_YC
 * @see fk_YC
 * @see Cache_YC
 */
double C11Y(double l, int maxIter, double rtol, double atol) {
    // Look for cached value
    auto& cache = Cache_YC::getInstance().CY_cache;
    auto key = std::make_pair(true, l);
    auto it = cache.find(key);
    if(it != cache.end()) {
        return it->second;
    }

    // Find or compute g-coefficients
    auto& g_cache = Cache_YC::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YC(l)).first;
    }
    const auto& g = g_it->second;

    double result = 0.0;
    double abs_sum = 0.0;

    double one_plus_l = 1.0 + l;
    double one_plus_l_inv = 1.0 / one_plus_l;

    for(int m = 2; m < 2*maxIter; m += 2) {

        double m_inv = 1.0/m;
        double m_plus_2_inv = 1.0/(m + 2.0);

        double term = (std::pow(0.5 * one_plus_l_inv, m) * fk_YC(m, l) - 
                        2.0 * m_inv * g.g2 + 4.0 * m_inv * m_plus_2_inv * g.g3);


        result += term;
        abs_sum += std::abs(term);

        if(m > 2 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: C11Y did not converge within " << maxIter << " iterations." << std::endl;
        }

    }

    result = 1 - g.g3 + result;

    // Store result in cache
    Cache_YC::getInstance().insertWithCheck(cache, key, result);

    return result;
}


double YC11_lubrication(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_YC::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YC(l)).first;
    }
    const auto& g = g_it->second;

    double ksi = s - 2.0;
    double ksi_inv = 1.0/ksi;
    double log_ksi_inv = std::log(ksi_inv);

    double result = g.g2 * log_ksi_inv + g.g3 * ksi * log_ksi_inv;
    result += C11Y(l, maxIter, rtol, atol);

    return result;
}

double YC11_nf(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_YC::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YC(l)).first;
    }
    const auto& g = g_it->second;
    
    double result = 0.0;
    double two_over_s = 2.0/s;
    double two_over_s_squared = two_over_s * two_over_s;
    double abs_sum = 0.0;
    double one_plus_l = 1.0 + l;
    double one_plus_l_inv = 1.0/one_plus_l;
    
    for(int m = 2; m < 2*maxIter; m += 2) {
        double m1_inv = (m == 2) ? -0.5 : 1.0/(m - 2.0);
        double m_inv = 1.0/m;

        double term = (std::pow(0.5 * one_plus_l_inv, m) * fk_YC(m, l) - 
                        2.0 * m_inv * g.g2 + 4.0 * m_inv * m1_inv * g.g3) * std::pow(two_over_s, m);
        
        result += term;
        abs_sum += std::abs(term);
        if(m > 2 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: YC11_nf did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    double arg = 1.0 - two_over_s_squared;
    double log_arg = std::log(arg);

    result = -g.g2 * log_arg - g.g3 * arg * log_arg + fk_YC(0, l) + result;

    return result;
}


double YC11_ff(double s, double l, int maxIter, double rtol, double atol) {
    double result = 0.0;
    double abs_sum = 0.0;
    double l_plus_1_inv = 1.0/(1.0 + l);
    double s_inv = 1.0/s;
    
    for(int k = 0; k < maxIter; k++) {
        double term = fk_YC(2*k, l) * std::pow(l_plus_1_inv * s_inv, 2 * k);
        result += term;
        abs_sum += std::abs(term);
        
        if(k > 0 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (k == maxIter - 1) {
            std::cerr << "Warning: YC11_ff did not converge within " << maxIter << " iterations." << std::endl;
        }

    }
    
    return result;
}

// ---------------------- YC12 functions ----------------------

/**
 * @brief Computes the C12Y coefficient for the two-body problem.
 * 
 * @details
 * This function calculates the C12Y coefficient using a series expansion with cached values
 * for optimization. It utilizes g-coefficients and implements a convergence check based on
 * relative and absolute tolerances. See Townsend (2023) for the formula.
 *
 * @param l The ratio of the two sphere radii
 * @param maxIter Maximum number of iterations for the series expansion
 * @param rtol Relative tolerance for convergence check
 * @param atol Absolute tolerance for convergence check
 * 
 * @return The computed C12Y coefficient value
 * 
 * @note The function uses caching mechanisms through Cache_YC for both C12Y values
 *       and g-coefficients to improve performance on repeated calls
 * @warning May print a warning to stderr if convergence is not reached within maxIter iterations
 * 
 * @see GCoeffs_YC
 * @see fk_YC
 * @see Cache_YC
 */
double C12Y(double l, int maxIter, double rtol, double atol) {
    // Look for cached value
    auto& cache = Cache_YC::getInstance().CY_cache;
    auto key = std::make_pair(false, l);
    auto it = cache.find(key);
    if(it != cache.end()) {
        return it->second;
    }

    // Find or compute g-coefficients
    auto& g_cache = Cache_YC::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YC(l)).first;
    }
    const auto& g = g_it->second;

    double result = 0.0;
    double abs_sum = 0.0;

    double one_plus_l = 1.0 + l;
    double one_plus_l_inv = 1.0 / one_plus_l;
    double one_plus_l_inv_cubed = one_plus_l_inv * one_plus_l_inv * one_plus_l_inv;

    for(int m = 1; m < 2*maxIter; m += 2) {
        double m_inv = 1.0/m;
        double m_plus_2_inv = 1.0/(m + 2.0);

        double term = (8.0 * one_plus_l_inv_cubed * std::pow(0.5 * one_plus_l_inv, m) * fk_YC(m, l) - 
                        2.0 * m_inv * g.g4 + 4.0 * m_inv * m_plus_2_inv * g.g5);

        result += term;
        abs_sum += std::abs(term);

        if(m > 1 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: C12Y did not converge within " << maxIter << " iterations." << std::endl;
        }

    }

    result = 2 * g.g4 * std::log(2.0) - 2 * g.g5 + result;

    // Store result in cache
    Cache_YC::getInstance().insertWithCheck(cache, key, result);

    return result;
}

double YC12_lubrication(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_YC::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YC(l)).first;
    }
    const auto& g = g_it->second;

    double ksi = s - 2.0;
    double ksi_inv = 1.0/ksi;
    double log_ksi_inv = std::log(ksi_inv);
    double one_plus_l = 1.0 + l;
    double one_plus_l_cube = one_plus_l * one_plus_l * one_plus_l;

    double result = g.g4 * log_ksi_inv + g.g5 * ksi * log_ksi_inv;
    result += C12Y(l, maxIter, rtol, atol);
    result *= one_plus_l_cube * 0.125;

    return result;
}


double YC12_nf(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_YC::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YC(l)).first;
    }
    const auto& g = g_it->second;
    
    double result = 0.0;
    double one_over_s = 1.0/s;
    double two_over_s = 2.0 * one_over_s;
    double two_over_s_squared = two_over_s * two_over_s;
    double one_plus_l = 1.0 + l;
    double one_plus_l_inv = 1.0/one_plus_l;
    double abs_sum = 0.0;

    
    for(int m = 1; m < 2*maxIter; m += 2) {
        double m1_inv = (m == 2) ? -0.5 : 1.0/(m - 2.0);
        double m_inv = 1.0/m;

        double term = (std::pow(0.5 * one_plus_l_inv, m) * fk_YC(m, l) - 
                        2.0 * m_inv * g.g4 + 4.0 * m_inv * m1_inv * g.g5) * std::pow(two_over_s, m);
        
        result += term;
        abs_sum += std::abs(term);
        if(m > 1 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: YC12_nf did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    double arg = 1.0 - two_over_s_squared;
    double log_term = std::log((s+2.0)/(s-2.0));

    result = g.g4 * log_term + g.g5 * arg * log_term + 4.0 * g.g5 * one_over_s + result;

    return result;
}

double YC12_ff(double s, double l, int maxIter, double rtol, double atol) {
    double result = 0.0;
    double abs_sum = 0.0;
    double l_plus_1_inv = 1.0/(1.0 + l);
    double s_inv = 1.0/s;
    
    for(int k = 0; k < maxIter; k++) {
        double term = fk_YC(2*k + 1, l) * std::pow(l_plus_1_inv * s_inv, 2.0*k + 1.0);
        result += term;
        abs_sum += std::abs(term);
        
        if(k > 0 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (k == maxIter - 1) {
            std::cerr << "Warning: YC12_ff did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    return result;
}


// ---------------------- Main functions ----------------------
double YC11(double s, double l, 
            double lubr_cutoff = DEFAULT_LUBR_CUTOFF,
            double cutoff = DEFAULT_CUTOFF,
            int maxIter = DEFAULT_MAX_ITER,
            double rtol = DEFAULT_RTOL,
            double atol = DEFAULT_ATOL) {
    
    if (s < lubr_cutoff) {
        return YC11_lubrication(s, l, maxIter, rtol, atol);
    } else if (s < cutoff) {
        return YC11_nf(s, l, maxIter, rtol, atol);
    } else {
        return YC11_ff(s, l, maxIter, rtol, atol);
    }
}

double YC12(double s, double l, 
            double lubr_cutoff = DEFAULT_LUBR_CUTOFF,
            double cutoff = DEFAULT_CUTOFF,
            int maxIter = DEFAULT_MAX_ITER,
            double rtol = DEFAULT_RTOL,
            double atol = DEFAULT_ATOL) {
    
    if (s < lubr_cutoff) {
        return YC12_lubrication(s, l, maxIter, rtol, atol);
    } else if (s < cutoff) {
        return YC12_nf(s, l, maxIter, rtol, atol);
    } else {
        return YC12_ff(s, l, maxIter, rtol, atol);
    }
}

// Add these implementations before the PYBIND11_MODULE:
namespace yc_utils {
    void clear_cache() {
        Cache_YC::getInstance().clear();
#ifdef ENABLE_PROFILING
        Profiler::getData().clear();
#endif
    }
}
