#include "XA.h"
#include "cache_hash.h"

#include "common.h"


/**
 * @struct GCoeffs_XA
 * @brief Calculates g1, g2 and g3 coefficients for the XA scalar resistance functions.
 *
 * These calculations follow the formulas (3.19a-c) in the paper by Jeffrey and Onishi (1984).
 * They are intermediate calculations used in the XA11 and XA12 functions.
 *
 * @param l Radius ratio parameter used in calculations
 *
 * Members:
 * @var g1 First coefficient = \lambda ^2/(1+\lambda )^3
 * @var g2 Second coefficient = \lambda (1 + 7\lambda + \lambda^2)/(1+\lambda)^3
 * @var g3 Third coefficient = (1/42) (1 + 18 \lambda - 29 \lambda^2 + 18 \lambda^3 + \lambda^4)/(1+\lambda)^3
 */
struct GCoeffs_XA {
    double g1, g2, g3;
    GCoeffs_XA(double l) {
        // Precompute common terms
        const double l2 = l * l;
        const double l3 = l2 * l;
        const double l4 = l2 * l2;
        const double one_plus_l = 1.0 + l;
        const double l_plus_1_cubed_inv = 1.0/(one_plus_l * one_plus_l * one_plus_l);
        
        // Calculate coefficients using precomputed terms
        g1 = 2.0 * l2 * l_plus_1_cubed_inv;
        g2 = 0.2 * l * (1.0 + 7.0*l + l2) * l_plus_1_cubed_inv;
        g3 = 0.02380952 * (1.0 + 18.0*l - 29.0*l2 + 18.0*l3 + l4) * l_plus_1_cubed_inv;
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
 * - AX_cache: Combined cache for A11X and A12X functions, indexed by type flag and l parameter
 * - g_coeffs_cache: Stores G-coefficients indexed by double
 * - comb_cache: Stores combinatorial values indexed by integer pairs
 * 
 * The singleton pattern ensures that only one instance of the cache exists throughout
 * the program's lifetime, providing a centralized storage for computed values.
 * 
 * Usage:
 * @code
 * Cache_XA& cache = Cache_XA::getInstance();
 * cache.clear(); // To clear all cached values
 * @endcode
 */
struct Cache_XA {
    static const size_t MAX_CACHE_SIZE = 100000000000; // Maximum entries per cache

    std::unordered_map<std::tuple<int,int,int>, double, tuple_int3_hash> P_cache; // (n, p, q) -> value
    std::unordered_map<std::tuple<int,int,int>, double, tuple_int3_hash> V_cache; // (n, p, q) -> value
    std::unordered_map<std::pair<int, double>, double, pair_int_double_hash> fk_cache; // (k, l) -> value
    std::unordered_map<std::pair<bool, double>, double, pair_bool_double_hash> AX_cache; // bool isXA11, double l
    std::unordered_map<double, GCoeffs_XA> g_coeffs_cache; // l -> GCoeffs_XA
    
    static Cache_XA& getInstance() {
        static Cache_XA instance;
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

        // fk_cache memory
        total_bytes += sizeof(std::pair<int, double>) * fk_cache.size() * 2;
        total_bytes += sizeof(double) * fk_cache.size();
        total_entries += fk_cache.size();

        // AX_cache memory
        total_bytes += sizeof(std::pair<bool, double>) * AX_cache.size() * 2;
        total_bytes += sizeof(double) * AX_cache.size();
        total_entries += AX_cache.size();

        // g_coeffs_cache memory
        total_bytes += sizeof(double) * g_coeffs_cache.size() * 4; // key + 3 coeffs
        total_entries += g_coeffs_cache.size();

        MemoryProfiler::record("Cache_XA", total_bytes, total_entries);
    }
#endif

    void clear() {
        P_cache.clear();
        V_cache.clear();
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
            // For unordered_map, remove half the elements
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
    Cache_XA() { 
#ifdef ENABLE_PROFILING
        updateMemoryStats(); 
#endif
    }
};



/**
 * @brief Calculates the pure P_{n,p,q} function for the XA scalar function in a two-body resistance problem.
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
 * @see V_npq_XA Function used in the recursive formula
 * @see Cache_XA::P_cache Cache storage for memoization
 */
double P_npq_XA(int n, int p, int q) {
#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    auto& cache = Cache_XA::getInstance().P_cache;
    auto key = std::make_tuple(n, p, q);
    auto it = cache.find(key);
    if(it != cache.end()) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("P_npq_XA", duration, cache_hit);
#endif
        return it->second;
    }

    if(p == 0 && q == 0) {
        Cache_XA::getInstance().insertWithCheck(cache, key, (n == 1) ? 1.0 : 0.0);
#ifdef ENABLE_PROFILING
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("P_npq_XA", duration, cache_hit);
#endif
        return cache[key];
    }

    double P = 0.0;
    double pre_term1 = n * (2.0*n + 1.0) / (2.0*(n + 1.0));
    double pre_term2 = -n * (2.0*n - 1.0) / (2.0*(n + 1.0));
    double pre_term3 = -n * (4.0*n*n - 1.0) / (2.0*(n + 1.0));
    for(int s = 1; s <= q; ++s) {
        double term1 = pre_term1 * (2.0*n*s - n - s + 2.0) / ((2.0*s - 1.0)*(n + s));
        double term2 = pre_term2;
        double term3 = pre_term3 / (2.0*s + 1.0);
        P += comb(n + s, n) * (
            term1 * P_npq_XA(s, q - s, p - n + 1) +
            term2 * P_npq_XA(s, q - s, p - n - 1) +
            term3 * V_npq_XA(s, q - s - 2, p - n + 1)
        );
    }

    Cache_XA::getInstance().insertWithCheck(cache, key, P);
#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("P_npq_XA", duration, cache_hit);
#endif
    return P;
}



/**
 * @brief Calculates the V_{n,p,q} pure function for the XA scalar function in a two-body resistance problem.
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
 * @note The function uses Cache_XA for memoization
 * @note Depends on P_npq_XA function and comb (combination) function
 * 
 * @see P_npq_XA(), comb()
 */
double V_npq_XA(int n, int p, int q) {
#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    auto& cache = Cache_XA::getInstance().V_cache;
    auto key = std::make_tuple(n, p, q);
    auto it = cache.find(key);
    if(it != cache.end()) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("V_npq_XA", duration, cache_hit);
#endif
        return it->second;
    }

    if(p == 0 && q == 0) {
        Cache_XA::getInstance().insertWithCheck(cache, key, (n == 1) ? 1.0 : 0.0);
#ifdef ENABLE_PROFILING
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("V_npq_XA", duration, cache_hit);
#endif
        return cache[key];
    }

    double term = 2.0 * n / ((n + 1.0)*(2.0*n + 3.0));
    double V = 0.0;
    for(int s = 1; s <= q; ++s) {
        V += comb(n + s, n) * P_npq_XA(s, q - s, p - n - 1);
    }

    V = P_npq_XA(n, p, q) - term * V;
    Cache_XA::getInstance().insertWithCheck(cache, key, V);
#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("V_npq_XA", duration, cache_hit);
#endif
    return V;
}



/**
 * @brief Calculates the fk function for the XA scalar function in a two-body resistance problem.
 * 
 * @details
 * This function implements a recursive calculation of fk with caching mechanism
 * for improved performance. It uses the following formula:
 * fk = 2^k * sum(P_1(k-q,q) * l^q) for q from 0 to k
 * 
 * @param k Integer parameter representing the order of the function
 * @param l Double parameter representing the ratio of sphere radii
 * 
 * @return Double value of the calculated fk function
 * 
 * @note The function uses a global cache to store previously calculated values
 * @note For k=0, only P_100 is needed
 * @note For k>0, it tries to use previously calculated k-1 terms if available
 * 
 * @see P_npq_XA, Cache_XA
 */
double fk_XA(int k, double l) {
#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    auto& cache = Cache_XA::getInstance().fk_cache;
    auto key = std::make_pair(k, l);
    auto it = cache.find(key);
    if(it != cache.end()) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("fk_XA", duration, cache_hit);
#endif
        return it->second;
    }

    double result = 0.0;
    double pow_2k = std::pow(2.0, k);
    

    for(int q = 0; q <= k; ++q) {
        // #ifdef _OPENMP
        // std::cout << "[fk_XA debug] Thread " << omp_get_thread_num() 
        //           << " handling q = " << q 
        //           << std::endl;
        // #endif
        result += P_npq_XA(1, k - q, q) * std::pow(l, q);
    }
    result *= pow_2k;

    Cache_XA::getInstance().insertWithCheck(cache, key, result);
#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("fk_XA", duration, cache_hit);
#endif
    return result;
}

// ---------------------- XA11 functions ----------------------

// Lubrication expression

/**
 * @brief Computes the A11X coefficient for the two-body problem.
 * 
 * @details
 * This function calculates the A11X coefficient using a series expansion with cached values
 * for optimization. It utilizes g-coefficients and implements a convergence check based on
 * relative and absolute tolerances. See equation 14a of Townsend (2023) for the formula.
 *
 * @param l The ratio of the two sphere radii
 * @param maxIter Maximum number of iterations for the series expansion
 * @param rtol Relative tolerance for convergence check
 * @param atol Absolute tolerance for convergence check
 * 
 * @return The computed A11X coefficient value
 * 
 * @note The function uses caching mechanisms through Cache_XA for both A11X values
 *       and g-coefficients to improve performance on repeated calls
 * @warning May print a warning to stderr if convergence is not reached within maxIter iterations
 * 
 * @see GCoeffs_XA
 * @see fk_XA
 * @see Cache_XA
 */
double A11X(double l, int maxIter, double rtol, double atol) {
    // Look for cached value
    auto& cache = Cache_XA::getInstance().AX_cache;
    auto key = std::make_pair(true, l); // Already using make_pair, no change needed
    auto it = cache.find(key);
    if(it != cache.end()) {
        return it->second;
    }

    // Find or compute g-coefficients
    auto& g_cache = Cache_XA::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_XA(l)).first;
    }
    const auto& g = g_it->second;

    double result = 0.0;
    double abs_sum = 0.0;

    for(int m = 2; m < 2*maxIter; m += 2) {
        double term = (std::pow(2.0, -m) * std::pow(1.0 + l, -m) * fk_XA(m, l) - 
                       g.g1 - 2.0/m * g.g2 + 4.0/(m * (m + 2.0)) * g.g3);
        result += term;
        abs_sum += std::abs(term);

        if(m > 2 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: A11X did not converge within " << maxIter << " iterations." << std::endl;
        }

    }

    result = 1 - 0.25 * g.g1 - g.g3 + result;

    // Store result in cache
    Cache_XA::getInstance().insertWithCheck(cache, key, result);

    return result;
}


double XA11_lubrication(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_XA::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_XA(l)).first;
    }
    const auto& g = g_it->second;

    double ksi = s - 2.0;
    double ksi_inv = 1.0/ksi;
    double log_ksi_inv = std::log(ksi_inv);

    double result = g.g1 * ksi_inv + g.g2 * log_ksi_inv + g.g3 * ksi * log_ksi_inv;
    result += A11X(l, maxIter, rtol, atol);

    return result;
}

double XA11_nf(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_XA::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_XA(l)).first;
    }
    const auto& g = g_it->second;
    
    double result = 0.0;
    double two_over_s = 2.0/s;
    double two_over_s_squared = two_over_s * two_over_s;
    double abs_sum = 0.0;
    
    for(int m = 2; m < 2*maxIter; m += 2) {
        double m1 = (m == 2) ? -2.0 : (m-2.0);
        double term = (std::pow(2.0, -m) * std::pow(1.0 + l, -m) * fk_XA(m, l) - 
                      g.g1 - 2.0/m * g.g2 + 4.0/(m * m1) * g.g3) * std::pow(two_over_s, m);
        
        result += term;
        abs_sum += std::abs(term);
        if(m > 2 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: XA11_nf did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    double arg = 1.0 - two_over_s_squared;

    result = g.g1/arg - g.g2*std::log(arg) - g.g3*arg*std::log(arg) + fk_XA(0, l) - g.g1 + result;

    return result;
}


double XA11_ff(double s, double l, int maxIter, double rtol, double atol) {
    double result = 0.0;
    double abs_sum = 0.0;
    double l_plus_1_inv = 1.0/(1.0 + l);
    double s_inv_squared = 1.0/(s*s);
    
    for(int k = 0; k < maxIter; k++) {
        double term = fk_XA(2*k, l) * std::pow(l_plus_1_inv, 2.0*k) * std::pow(s_inv_squared, k);
        result += term;
        abs_sum += std::abs(term);
        
        if(k > 0 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (k == maxIter - 1) {
            std::cerr << "Warning: XA11_ff did not converge within " << maxIter << " iterations." << std::endl;
        }

    }
    
    return result;
}

// ---------------------- XA12 functions ----------------------

/**
 * @brief Computes the A12X coefficient for the two-body problem.
 * 
 * @details
 * This function calculates the A12X coefficient using a series expansion with cached values
 * for optimization. It utilizes g-coefficients and implements a convergence check based on
 * relative and absolute tolerances. See equation 15 of Townsend (2023) for the formula.
 *
 * @param l The ratio of the two sphere radii
 * @param maxIter Maximum number of iterations for the series expansion
 * @param rtol Relative tolerance for convergence check
 * @param atol Absolute tolerance for convergence check
 * 
 * @return The computed A12X coefficient value
 * 
 * @note The function uses caching mechanisms through Cache_XA for both A12X values
 *       and g-coefficients to improve performance on repeated calls
 * @warning May print a warning to stderr if convergence is not reached within maxIter iterations
 * 
 * @see GCoeffs_XA
 * @see fk_XA
 * @see Cache_XA
 */
double A12X(double l, int maxIter, double rtol, double atol) {
    // Look for cached value
    auto& cache = Cache_XA::getInstance().AX_cache;
    auto key = std::make_pair(false, l); // Already using make_pair, no change needed
    auto it = cache.find(key);
    if(it != cache.end()) {
        return it->second;
    }

    // Find or compute g-coefficients
    auto& g_cache = Cache_XA::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_XA(l)).first;
    }
    const auto& g = g_it->second;

    double result = 0.0;
    double abs_sum = 0.0;

    for(int m = 1; m < 2*maxIter; m += 2) {
        double term = (std::pow(2.0, -m) * std::pow(1.0 + l, -m) * fk_XA(m, l) - 
                       g.g1 - 2.0/m * g.g2 + 4.0/(m * (m+2.0)) * g.g3);

        result += term;
        abs_sum += std::abs(term);

        if(m > 1 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: A12X did not converge within " << maxIter << " iterations." << std::endl;
        }

    }

    result = 0.25 * g.g1 + 2.0 * g.g2 * std::log(2.0) - 2.0 * g.g3 + result;

    // We ommit the -1/2 (1 + l) term in the original formula because it's included in the XA12_lubrication function

    // Store result in cache
    Cache_XA::getInstance().insertWithCheck(cache, key, result);

    return result;
}

double XA12_lubrication(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_XA::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_XA(l)).first;
    }
    const auto& g = g_it->second;

    double ksi = s - 2.0;
    double ksi_inv = 1.0/ksi;
    double log_ksi_inv = std::log(ksi_inv);

    double result = g.g1 * ksi_inv + g.g2 * log_ksi_inv + g.g3 * ksi * log_ksi_inv;
    result += A12X(l, maxIter, rtol, atol);

    return -result;
}


double XA12_nf(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_XA::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_XA(l)).first;
    }
    const auto& g = g_it->second;
    
    double result = 0.0;
    double two_over_s = 2.0/s;
    double two_over_s_squared = two_over_s * two_over_s;
    double abs_sum = 0.0;
    
    for(int m = 1; m < 2*maxIter; m += 2) {
        double m1 = (m == 2) ? -2.0 : (m-2.0);
        double term = (std::pow(2.0, -m) * std::pow(1.0 + l, -m) * fk_XA(m, l) - 
                      g.g1 - 2.0/m * g.g2 + 4.0/(m * m1) * g.g3) * std::pow(two_over_s, m);
        
        result += term;
        abs_sum += std::abs(term);
        if(m > 1 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: XA12_nf did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    double arg = 1.0 - two_over_s_squared;
    double log_term = std::log((s+2.0)/(s-2.0));
    result = two_over_s * g.g1/arg + g.g2*log_term + g.g3*arg*log_term + 4.0*g.g3/s + result;
    
    return -result;
}

double XA12_ff(double s, double l, int maxIter, double rtol, double atol) {
    double result = 0.0;
    double abs_sum = 0.0;
    double l_plus_1_inv = 1.0/(1.0 + l);
    double s_inv = 1.0/s;
    
    for(int k = 0; k < maxIter; k++) {
        double term = fk_XA(2*k+1, l) * std::pow(l_plus_1_inv, 2.0*k+1.0) * std::pow(s_inv, 2.0*k+1.0);
        result += term;
        abs_sum += std::abs(term);
        
        if(k > 0 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (k == maxIter - 1) {
            std::cerr << "Warning: XA12_ff did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    return -result;
}


// ---------------------- Main functions ----------------------
double XA11(double s, double l, 
            double lubr_cutoff = DEFAULT_LUBR_CUTOFF,
            double cutoff = DEFAULT_CUTOFF,
            int maxIter = DEFAULT_MAX_ITER,
            double rtol = DEFAULT_RTOL,
            double atol = DEFAULT_ATOL) {
    
    if (s < lubr_cutoff) {
        return XA11_lubrication(s, l, maxIter, rtol, atol);
    } else if (s < cutoff) {
        return XA11_nf(s, l, maxIter, rtol, atol);
    } else {
        return XA11_ff(s, l, maxIter, rtol, atol);
    }
}

double XA12(double s, double l, 
            double lubr_cutoff = DEFAULT_LUBR_CUTOFF,
            double cutoff = DEFAULT_CUTOFF,
            int maxIter = DEFAULT_MAX_ITER,
            double rtol = DEFAULT_RTOL,
            double atol = DEFAULT_ATOL) {
    
    if (s < lubr_cutoff) {
        return XA12_lubrication(s, l, maxIter, rtol, atol);
    } else if (s < cutoff) {
        return XA12_nf(s, l, maxIter, rtol, atol);
    } else {
        return XA12_ff(s, l, maxIter, rtol, atol);
    }
}

// Add these implementations before the PYBIND11_MODULE:
namespace xa_utils {
    void clear_cache() {
        Cache_XA::getInstance().clear();
#ifdef ENABLE_PROFILING
        Profiler::getData().clear();
#endif
    }
}