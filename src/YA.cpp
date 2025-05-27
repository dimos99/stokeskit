#include "YA.h"
#include "cache_hash.h"

#include "common.h"


/**
 * @struct GCoeffs_YA
 * @brief Calculates g1, g2 and g3 coefficients for the YA scalar resistance functions.
 *
 * These calculations follow the formulas in the paper by Jeffrey and Onishi (1984).
 * They are intermediate calculations used in the YA11 and YA12 functions.
 *
 * @param l Radius ratio parameter used in calculations
 *
 * Members:
 * @var g2 Second coefficient
 * @var g3 Third coefficient
 */
struct GCoeffs_YA {
    double g2, g3;
    GCoeffs_YA(double l) {
        // Precompute common terms
        const double l2 = l * l;
        const double l3 = l2 * l;
        const double l4 = l2 * l2;
        const double one_plus_l = 1.0 + l;
        const double l_plus_1_cubed_inv = 1.0/(one_plus_l * one_plus_l * one_plus_l);
        
        // Calculate coefficients using precomputed terms
        g2 = 0.26666667 * l * (2.0 + l + 2*l2) * l_plus_1_cubed_inv;
        g3 = 0.00533333 * (16.0 - 45.0 * l + 58.0 * l2 - 45.0 * l3 + 16.0 * l4) * l_plus_1_cubed_inv;
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
 * - AX_cache: Combined cache for A11Y and A12Y functions, indexed by type flag and l parameter
 * - g_coeffs_cache: Stores G-coefficients indexed by double
 * - comb_cache: Stores combinatorial values indexed by integer pairs
 * 
 * The singleton pattern ensures that only one instance of the cache exists throughout
 * the program's lifetime, providing a centralized storage for computed values.
 * 
 * Usage:
 * @code
 * Cache_YA& cache = Cache_YA::getInstance();
 * cache.clear(); // To clear all cached values
 * @endcode
 */
struct Cache_YA {
    static const size_t MAX_CACHE_SIZE = 100000000000; // Updated size to match XA

    std::unordered_map<std::tuple<int,int,int>, double, tuple_int3_hash> P_cache; // (n, p, q) -> value
    std::unordered_map<std::tuple<int,int,int>, double, tuple_int3_hash> V_cache; // (n, p, q) -> value
    std::unordered_map<std::tuple<int,int,int>, double, tuple_int3_hash> Q_cache; // (n, p, q) -> value
    std::unordered_map<std::pair<int, double>, double, pair_int_double_hash> fk_cache; // (k, l) -> value
    std::unordered_map<std::pair<bool, double>, double, pair_bool_double_hash> AX_cache; // bool isYA11, double l
    std::unordered_map<double, GCoeffs_YA> g_coeffs_cache; // l -> GCoeffs_YA
    
    static Cache_YA& getInstance() {
        static Cache_YA instance;
        return instance;
    }

#ifdef ENABLE_PROFILING
    void updateMemoryStats() {
        size_t total_bytes = 0;
        size_t total_entries = 0;

        // P_cache memory
        total_bytes += sizeof(std::tuple<int,int,int>) * P_cache.size() * 2;
        total_bytes += sizeof(double) * P_cache.size();
        total_entries += P_cache.size();

        // V_cache memory
        total_bytes += sizeof(std::tuple<int,int,int>) * V_cache.size() * 2;
        total_bytes += sizeof(double) * V_cache.size();
        total_entries += V_cache.size();

        // Q_cache memory
        total_bytes += sizeof(std::tuple<int,int,int>) * Q_cache.size() * 2;
        total_bytes += sizeof(double) * Q_cache.size();
        total_entries += Q_cache.size();

        // fk_cache memory
        total_bytes += sizeof(std::pair<int, double>) * fk_cache.size() * 2;
        total_bytes += sizeof(double) * fk_cache.size();
        total_entries += fk_cache.size();

        // AX_cache memory
        total_bytes += sizeof(std::pair<bool, double>) * AX_cache.size() * 2;
        total_bytes += sizeof(double) * AX_cache.size();
        total_entries += AX_cache.size();

        // g_coeffs_cache memory
        total_bytes += sizeof(double) * g_coeffs_cache.size() * 3; // key + 2 coeffs
        total_entries += g_coeffs_cache.size();

        MemoryProfiler::record("Cache_YA", total_bytes, total_entries);
    }
#endif

    void clear() {
        P_cache.clear();
        V_cache.clear();
        Q_cache.clear();
        fk_cache.clear();
        g_coeffs_cache.clear();
#ifdef ENABLE_PROFILING
        updateMemoryStats();
#endif
    }

    template<typename MapType>
    void checkAndTrimCache(MapType& cache) {
        if (cache.size() > MAX_CACHE_SIZE) {
            std::cerr << "Warning: Cache size exceeded maximum limit. Clearing cache..." << std::endl;
            auto it = cache.begin();
            size_t to_remove = cache.size() / 2;
            for (size_t i = 0; i < to_remove && it != cache.end(); ++i) {
                it = cache.erase(it);
            }
        }
    }

    template<typename K, typename V, typename H>
    void insertWithCheck(std::unordered_map<K,V,H>& cache, const K& key, const V& value) {
        cache[key] = value;
        checkAndTrimCache(cache);
#ifdef ENABLE_PROFILING
        updateMemoryStats();
#endif
    }

private:
    Cache_YA() {
#ifdef ENABLE_PROFILING
        updateMemoryStats(); 
#endif
    }
};


/**
 * @brief Calculates the pure P_{n,p,q} function for the YA scalar function in a two-body resistance problem.
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
 * @see V_npq_YA Function used in the recursive formula
 * @see Q_npq_YA Function used in the recursive formula
 * @see Cache_YA::P_cache Cache storage for memoization
 */
double P_npq_YA(int n, int p, int q) {
#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    auto& cache = Cache_YA::getInstance().P_cache;
    auto key = std::make_tuple(n, p, q);
    auto it = cache.find(key);
    if(it != cache.end()) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("P_npq_YA", duration, cache_hit);
#endif
        return it->second;
    }

    // Base case
    if(p == 0 && q == 0) {
        Cache_YA::getInstance().insertWithCheck(cache, key, (n == 1) ? 1.0 : 0.0);
#ifdef ENABLE_PROFILING
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("P_npq_YA", duration, cache_hit);
#endif
        return cache[key];
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
            term1 * P_npq_YA(s, q - s, p - n + 1) +
            term2 * P_npq_YA(s, q - s, p - n - 1) +
            term3 * V_npq_YA(s, q - s - 2, p - n + 1) +
            term4 * Q_npq_YA(s, q - s - 1, p - n + 1)
        );
    }

    Cache_YA::getInstance().insertWithCheck(cache, key, P);
#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("P_npq_YA", duration, cache_hit);
#endif
    return P;
}



/**
 * @brief Calculates the V_{n,p,q} pure function for the YA scalar function in a two-body resistance problem.
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
 * @note The function uses Cache_YA for memoization
 * @note Depends on P_npq_YA function and comb (combination) function
 * @note This function is wrong in Jeffrey & Onishi (1984) paper, the correct formula is in Townsend (2023)
 * 
 * @see P_npq_YA(), comb()
 */
double V_npq_YA(int n, int p, int q) {
#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    auto& cache = Cache_YA::getInstance().V_cache;
    auto key = std::make_tuple(n, p, q);
    auto it = cache.find(key);
    if(it != cache.end()) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("V_npq_YA", duration, cache_hit);
#endif
        return it->second;
    }

    if(p == 0 && q == 0) {
        Cache_YA::getInstance().insertWithCheck(cache, key, (n == 1) ? 1.0 : 0.0);
#ifdef ENABLE_PROFILING
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("V_npq_YA", duration, cache_hit);
#endif
        return cache[key];
    }

    double term = 2.0 * n / ((n + 1.0)*(2.0*n + 3.0));
    double V = 0.0;
    for(int s = 1; s <= q; ++s) {
        V += comb(n + s, n + 1) * P_npq_YA(s, q - s, p - n - 1);
    }

    V = P_npq_YA(n, p, q) + term * V;
    Cache_YA::getInstance().insertWithCheck(cache, key, V);
#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("V_npq_YA", duration, cache_hit);
#endif
    return V;
}


/** 
 * @brief Calculates the Q_{n,p,q} pure function for the YA scalar function in a two-body resistance problem.
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
 * @note The function uses Cache_YA for memoization
 * @note Depends on P_npq_YA function and comb (combination) function
 * 
 * @see P_npq_YA(), comb()
 */
double Q_npq_YA(int n, int p, int q) {
#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    auto& cache = Cache_YA::getInstance().Q_cache;
    auto key = std::make_tuple(n, p, q);
    auto it = cache.find(key);
    if(it != cache.end()) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("Q_npq_YA", duration, cache_hit);
#endif
        return it->second;
    }

    if(p == 0 && q == 0) {
        Cache_YA::getInstance().insertWithCheck(cache, key, 0.0);
#ifdef ENABLE_PROFILING
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("Q_npq_YA", duration, cache_hit);
#endif
        return cache[key];
    }

    double pre_term1 = 1.0 / (n + 1.0);
    double pre_term2 = - 3.0 / (2.0 * n * (n + 1.0));
    double Q = 0.0;
    for(int s = 1; s <= q; ++s) {
        Q += comb(n + s, n + 1) * (
            pre_term1 * s * Q_npq_YA(s, q - s - 1, p - n) +
            pre_term2 / s * P_npq_YA(s, q - s, p - n)
        );
    }

    Cache_YA::getInstance().insertWithCheck(cache, key, Q);

#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("Q_npq_YA", duration, cache_hit);
#endif
    return Q;
}


/**
 * @brief Calculates the fk function for the YA scalar function in a two-body resistance problem.
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
 * @see P_npq_YA, Cache_YA
 */
double fk_YA(int k, double l) {
#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    auto& cache = Cache_YA::getInstance().fk_cache;
    auto key = std::make_pair(k, l);
    auto it = cache.find(key);
    if(it != cache.end()) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("fk_YA", duration, cache_hit);
#endif
        return it->second;
    }

    double result = 0.0;
    double pow_2k = std::pow(2.0, k);
    

    for(int q = 0; q <= k; ++q) {
        // #ifdef _OPENMP
        // std::cout << "[fk_YA debug] Thread " << omp_get_thread_num() 
        //           << " handling q = " << q 
        //           << std::endl;
        // #endif
        result += P_npq_YA(1, k - q, q) * std::pow(l, q);
    }
    result *= pow_2k;

    Cache_YA::getInstance().insertWithCheck(cache, key, result);
#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("fk_YA", duration, cache_hit);
#endif
    return result;
}

// ---------------------- YA11 functions ----------------------

// Lubrication expression

/**
 * @brief Computes the A11Y coefficient for the two-body problem.
 * 
 * @details
 * This function calculates the A11Y coefficient using a series expansion with cached values
 * for optimization. It utilizes g-coefficients and implements a convergence check based on
 * relative and absolute tolerances. See equation 18a of Townsend (2023) for the formula.
 *
 * @param l The ratio of the two sphere radii
 * @param maxIter Maximum number of iterations for the series expansion
 * @param rtol Relative tolerance for convergence check
 * @param atol Absolute tolerance for convergence check
 * 
 * @return The computed A11Y coefficient value
 * 
 * @note The function uses caching mechanisms through Cache_YA for both A11Y values
 *       and g-coefficients to improve performance on repeated calls
 * @warning May print a warning to stderr if convergence is not reached within maxIter iterations
 * 
 * @see GCoeffs_YA
 * @see fk_YA
 * @see Cache_YA
 */
double A11Y(double l, int maxIter, double rtol, double atol) {
    // Look for cached value
    auto& cache = Cache_YA::getInstance().AX_cache;
    auto key = std::make_pair(true, l);
    auto it = cache.find(key);
    if(it != cache.end()) {
        return it->second;
    }

    // Find or compute g-coefficients
    auto& g_cache = Cache_YA::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YA(l)).first;
    }
    const auto& g = g_it->second;

    double result = 0.0;
    double abs_sum = 0.0;

    for(int m = 2; m < 2*maxIter; m += 2) {
        double term = (std::pow(2.0 * (1.0 + l), -m) * fk_YA(m, l) - 
                       2.0/m * g.g2 + 4.0/(m * (m + 2.0)) * g.g3);
        result += term;
        abs_sum += std::abs(term);

        if(m > 2 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: A11Y did not converge within " << maxIter << " iterations." << std::endl;
        }

    }

    result = 1 - g.g3 + result;

    // Store result in cache
    Cache_YA::getInstance().insertWithCheck(cache, key, result);

    return result;
}


double YA11_lubrication(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_YA::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YA(l)).first;
    }
    const auto& g = g_it->second;

    double ksi = s - 2.0;
    double ksi_inv = 1.0/ksi;
    double log_ksi_inv = std::log(ksi_inv);

    double result = g.g2 * log_ksi_inv + g.g3 * ksi * log_ksi_inv;
    result += A11Y(l, maxIter, rtol, atol);

    return result;
}

double YA11_nf(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_YA::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YA(l)).first;
    }
    const auto& g = g_it->second;
    
    double result = 0.0;
    double two_over_s = 2.0/s;
    double two_over_s_squared = two_over_s * two_over_s;
    double abs_sum = 0.0;
    
    for(int m = 2; m < 2*maxIter; m += 2) {
        double m1 = (m == 2) ? -2.0 : (m-2.0);
        double term = (std::pow(2.0 * (1.0 + l), -m) * fk_YA(m, l) - 
                      2.0/m * g.g2 + 4.0/(m * m1) * g.g3) * std::pow(two_over_s, m);
        
        result += term;
        abs_sum += std::abs(term);
        if(m > 2 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: YA11_nf did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    double arg = 1.0 - two_over_s_squared;

    result = - (g.g2 + g.g3 * arg) * std::log(arg) + fk_YA(0, l) + result;

    return result;
}


double YA11_ff(double s, double l, int maxIter, double rtol, double atol) {
    double result = 0.0;
    double abs_sum = 0.0;
    double l_plus_1_inv = 1.0/(1.0 + l);
    double s_inv = 1.0/s;
    
    for(int k = 0; k < maxIter; k++) {
        double term = fk_YA(2*k, l) * std::pow(l_plus_1_inv * s_inv, 2 * k);
        result += term;
        abs_sum += std::abs(term);
        
        if(k > 0 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (k == maxIter - 1) {
            std::cerr << "Warning: YA11_ff did not converge within " << maxIter << " iterations." << std::endl;
        }

    }
    
    return result;
}

// ---------------------- YA12 functions ----------------------

/**
 * @brief Computes the A12Y coefficient for the two-body problem.
 * 
 * @details
 * This function calculates the A12Y coefficient using a series expansion with cached values
 * for optimization. It utilizes g-coefficients and implements a convergence check based on
 * relative and absolute tolerances. See Townsend (2023) for the formula.
 *
 * @param l The ratio of the two sphere radii
 * @param maxIter Maximum number of iterations for the series expansion
 * @param rtol Relative tolerance for convergence check
 * @param atol Absolute tolerance for convergence check
 * 
 * @return The computed A12Y coefficient value
 * 
 * @note The function uses caching mechanisms through Cache_YA for both A12Y values
 *       and g-coefficients to improve performance on repeated calls
 * @warning May print a warning to stderr if convergence is not reached within maxIter iterations
 * 
 * @see GCoeffs_YA
 * @see fk_YA
 * @see Cache_YA
 */
double A12Y(double l, int maxIter, double rtol, double atol) {
    // Look for cached value
    auto& cache = Cache_YA::getInstance().AX_cache;
    auto key = std::make_pair(false, l);
    auto it = cache.find(key);
    if(it != cache.end()) {
        return it->second;
    }

    // Find or compute g-coefficients
    auto& g_cache = Cache_YA::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YA(l)).first;
    }
    const auto& g = g_it->second;

    double result = 0.0;
    double abs_sum = 0.0;

    for(int m = 1; m < 2*maxIter; m += 2) {
        double term = (std::pow(2.0 * (1.0 + l), -m) * fk_YA(m, l) - 
                          2.0/m * g.g2 + 4.0/(m * (m + 2.0)) * g.g3);

        result += term;
        abs_sum += std::abs(term);

        if(m > 1 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: A12Y did not converge within " << maxIter << " iterations." << std::endl;
        }

    }

    result = 2 * g.g2 * std::log(2.0) - 2.0 * g.g3 + result;

    // We ommit the -1/2 (1 + l) term in the original formula because it's included in the YA12_lubrication function

    // Store result in cache
    Cache_YA::getInstance().insertWithCheck(cache, key, result);

    return result;
}

double YA12_lubrication(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_YA::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YA(l)).first;
    }
    const auto& g = g_it->second;

    double ksi = s - 2.0;
    double ksi_inv = 1.0/ksi;
    double log_ksi_inv = std::log(ksi_inv);

    double result = g.g2 * log_ksi_inv + g.g3 * ksi * log_ksi_inv;
    result += A12Y(l, maxIter, rtol, atol);

    return -result;
}


double YA12_nf(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_YA::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YA(l)).first;
    }
    const auto& g = g_it->second;
    
    double result = 0.0;
    double two_over_s = 2.0/s;
    double two_over_s_squared = two_over_s * two_over_s;
    double abs_sum = 0.0;
    
    for(int m = 1; m < 2*maxIter; m += 2) {
        double m1 = (m == 2) ? -2.0 : (m-2.0);
        double term = (std::pow(2.0 * (1.0 + l), -m) * fk_YA(m, l) - 
                      2.0/m * g.g2 + 4.0/(m * m1) * g.g3) * std::pow(two_over_s, m);
        
        result += term;
        abs_sum += std::abs(term);
        if(m > 1 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: YA12_nf did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    double arg = 1.0 - two_over_s_squared;
    double log_term = std::log((s+2.0)/(s-2.0));
    result = g.g2 * log_term + g.g3 *  arg * log_term + 2 * g.g3 * two_over_s + result;
    
    return -result;
}

double YA12_ff(double s, double l, int maxIter, double rtol, double atol) {
    double result = 0.0;
    double abs_sum = 0.0;
    double l_plus_1_inv = 1.0/(1.0 + l);
    double s_inv = 1.0/s;
    
    for(int k = 0; k < maxIter; k++) {
        double term = fk_YA(2*k+1, l) * std::pow(l_plus_1_inv * s_inv, 2.0*k+1.0);
        result += term;
        abs_sum += std::abs(term);
        
        if(k > 0 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (k == maxIter - 1) {
            std::cerr << "Warning: YA12_ff did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    return -result;
}


// ---------------------- Main functions ----------------------
double YA11(double s, double l, 
            double lubr_cutoff = DEFAULT_LUBR_CUTOFF,
            double cutoff = DEFAULT_CUTOFF,
            int maxIter = DEFAULT_MAX_ITER,
            double rtol = DEFAULT_RTOL,
            double atol = DEFAULT_ATOL) {
    
    if (s < lubr_cutoff) {
        return YA11_lubrication(s, l, maxIter, rtol, atol);
    } else if (s < cutoff) {
        return YA11_nf(s, l, maxIter, rtol, atol);
    } else {
        return YA11_ff(s, l, maxIter, rtol, atol);
    }
}

double YA12(double s, double l, 
            double lubr_cutoff = DEFAULT_LUBR_CUTOFF,
            double cutoff = DEFAULT_CUTOFF,
            int maxIter = DEFAULT_MAX_ITER,
            double rtol = DEFAULT_RTOL,
            double atol = DEFAULT_ATOL) {
    
    if (s < lubr_cutoff) {
        return YA12_lubrication(s, l, maxIter, rtol, atol);
    } else if (s < cutoff) {
        return YA12_nf(s, l, maxIter, rtol, atol);
    } else {
        return YA12_ff(s, l, maxIter, rtol, atol);
    }
}

// Add these implementations before the PYBIND11_MODULE:
namespace ya_utils {
    void clear_cache() {
        Cache_YA::getInstance().clear();
#ifdef ENABLE_PROFILING
        Profiler::getData().clear();
#endif
    }
}
