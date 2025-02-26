#include "YB.h"

#include "common.h"


/**
 * @struct GCoeffs_YB
 * @brief Calculates g1, g2 and g3 coefficients for the YB scalar resistance functions.
 *
 * These calculations follow the formulas in the paper by Jeffrey and Onishi (1984).
 * They are intermediate calculations used in the YB11 and YB12 functions.
 *
 * @param l Radius ratio parameter used in calculations
 *
 * Members:
 * @var g2 Second coefficient
 * @var g3 Third coefficient
 */
struct GCoeffs_YB {
    double g2, g3;
    GCoeffs_YB(double l) {
        // Precompute common terms
        const double l2 = l * l;
        const double l3 = l2 * l;
        const double one_plus_l = 1.0 + l;
        const double l_plus_1_square_inv = 1.0/(one_plus_l * one_plus_l);
        
        // Calculate coefficients using precomputed terms
        g2 = -0.20 * l * (4.0 + l) * l_plus_1_square_inv;
        g3 = -0.004 * (32.0 - 32.0 * l + 83.0 * l2 + 43.0 * l3) * l_plus_1_square_inv;
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
 * - BX_cache: Combined cache for B11Y and B12Y functions, indexed by type flag and l parameter
 * - g_coeffs_cache: Stores G-coefficients indexed by double
 * - comb_cache: Stores combinatorial values indexed by integer pairs
 * 
 * The singleton pattern ensures that only one instance of the cache exists throughout
 * the program's lifetime, providing a centralized storage for computed values.
 * 
 * Usage:
 * @code
 * Cache_YB& cache = Cache_YB::getInstance();
 * cache.clear(); // To clear all cached values
 * @endcode
 */
struct Cache_YB {
    static const size_t MAX_CACHE_SIZE = 100000000000; // Updated size to match XA

    std::map<std::pair<int, double>, double> fk_cache; // (k, l) -> value
    std::map<std::pair<bool, double>, double> BX_cache; // bool isYB11, double l
    std::map<double, GCoeffs_YB> g_coeffs_cache; // l -> GCoeffs_YB
    
    static Cache_YB& getInstance() {
        static Cache_YB instance;
        return instance;
    }
#ifdef ENABLE_PROFILING
    void updateMemoryStats() {
        size_t total_bytes = 0;
        size_t total_entries = 0;

        // fk_cache memory
        total_bytes += sizeof(std::pair<int, double>) * fk_cache.size() * 2;
        total_bytes += sizeof(double) * fk_cache.size();
        total_entries += fk_cache.size();

        // BX_cache memory
        total_bytes += sizeof(std::pair<bool, double>) * BX_cache.size() * 2;
        total_bytes += sizeof(double) * BX_cache.size();
        total_entries += BX_cache.size();

        // g_coeffs_cache memory
        total_bytes += sizeof(double) * g_coeffs_cache.size() * 3; // key + 2 coeffs
        total_entries += g_coeffs_cache.size();

        MemoryProfiler::record("Cache_YB", total_bytes, total_entries);
    }
#endif

    void clear() {
        fk_cache.clear();
        g_coeffs_cache.clear();
        BX_cache.clear();
#ifdef ENABLE_PROFILING        
        updateMemoryStats();
#endif
    }

    template<typename MapType>
    void checkAndTrimCache(MapType& cache) {
        if (cache.size() > MAX_CACHE_SIZE) {
            std::cerr << "Warning: Cache size exceeded maximum limit. Clearing cache..." << std::endl;
            auto it = cache.begin();
            std::advance(it, cache.size() / 2);
            cache.erase(it, cache.end());
        }
    }

    template<typename K, typename V>
    void insertWithCheck(std::map<K,V>& cache, const K& key, const V& value) {
        cache[key] = value;
        checkAndTrimCache(cache);
#ifdef ENABLE_PROFILING
        updateMemoryStats();
#endif
    }

private:
    Cache_YB() { 
#ifdef ENABLE_PROFILING
        updateMemoryStats();
#endif
    }
};


/**
 * @brief Calculates the pure P_{n,p,q} function for the YB scalar function in a two-body resistance problem.
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
 * @see V_npq_YB Function used in the recursive formula
 * @see Q_npq_YB Function used in the recursive formula
 * @see Cache_YB::P_cache Cache storage for memoization
 */
double P_npq_YB(int n, int p, int q) {
    return P_npq_YA(n, p, q);
}



/**
 * @brief Calculates the V_{n,p,q} pure function for the YB scalar function in a two-body resistance problem.
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
 * @note The function uses Cache_YB for memoization
 * @note Depends on P_npq_YB function and comb (combination) function
 * 
 * @see P_npq_YB(), comb()
 */
double V_npq_YB(int n, int p, int q) {
    return V_npq_YA(n, p, q);
}


/** 
 * @brief Calculates the Q_{n,p,q} pure function for the YB scalar function in a two-body resistance problem.
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
 * @note The function uses Cache_YB for memoization
 * @note Depends on P_npq_YB function and comb (combination) function
 * 
 * @see P_npq_YB(), comb()
 */
double Q_npq_YB(int n, int p, int q) {
    return Q_npq_YA(n, p, q);
}


/**
 * @brief Calculates the fk function for the YB scalar function in a two-body resistance problem.
 * 
 * @details
 * This function implements a recursive calculation of fk with caching mechanism
 * for improved performance. See Jeffrey and Onishi (1984).
 * 
 * @param k Integer parameter representing the order of the function
 * @param l Double parameter representing the ratio of sphere radii
 * 
 * @return Double value of the calculated fk function
 * 
 * @note The function uses a global cache to store previously calculated values
 * 
 * @see Q_npq_YB, Cache_YB
 */
double fk_YB(int k, double l) {
#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    auto& cache = Cache_YB::getInstance().fk_cache;
    auto key = std::make_pair(k, l);
    auto it = cache.find(key);
    if(it != cache.end()) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("fk_YB", duration, cache_hit);
#endif
        return it->second;
    }

    double result = 0.0;
    

    for(int q = 0; q <= k; ++q) {
        // #ifdef _OPENMP
        // std::cout << "[fk_YB debug] Thread " << omp_get_thread_num() 
        //           << " handling q = " << q 
        //           << std::endl;
        // #endif
        result += Q_npq_YB(1, k - q, q) * std::pow(l, q);
    }
    result *= std::pow(2.0, k+1);

    Cache_YB::getInstance().insertWithCheck(cache, key, result);
#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("fk_YB", duration, cache_hit);
#endif
    return result;
}

// ---------------------- YB11 functions ----------------------

// Lubrication expression

/**
 * @brief Computes the B11Y coefficient for the two-body problem.
 * 
 * @details
 * This function calculates the B11Y coefficient using a series expansion with cached values
 * for optimization. It utilizes g-coefficients and implements a convergence check based on
 * relative and absolute tolerances. See equation 18a of Townsend (2023) for the formula.
 *
 * @param l The ratio of the two sphere radii
 * @param maxIter Maximum number of iterations for the series expansion
 * @param rtol Relative tolerance for convergence check
 * @param atol Absolute tolerance for convergence check
 * 
 * @return The computed B11Y coefficient value
 * 
 * @note The function uses caching mechanisms through Cache_YB for both B11Y values
 *       and g-coefficients to improve performance on repeated calls
 * @warning May print a warning to stderr if convergence is not reached within maxIter iterations
 * 
 * @see GCoeffs_YB
 * @see fk_YB
 * @see Cache_YB
 */
double B11Y(double l, int maxIter, double rtol, double atol) {
    // Look for cached value
    auto& cache = Cache_YB::getInstance().BX_cache;
    auto key = std::make_pair(true, l);
    auto it = cache.find(key);
    if(it != cache.end()) {
        return it->second;
    }

    // Find or compute g-coefficients
    auto& g_cache = Cache_YB::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YB(l)).first;
    }
    const auto& g = g_it->second;

    double result = 0.0;
    double abs_sum = 0.0;

    for(int m = 1; m < 2*maxIter; m += 2) {
        double term = (std::pow(2.0 * (1.0 + l), -m) * fk_YB(m, l) - 
                       2.0/m * g.g2 + 4.0/(m * (m + 2.0)) * g.g3);
        result += term;
        abs_sum += std::abs(term);

        if(m > 2 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: B11Y did not converge within " << maxIter << " iterations." << std::endl;
        }

    }

    result = 2.0 * g.g2 * std::log(2.0) - 2.0 * g.g3 + result;

    // Store result in cache
    Cache_YB::getInstance().insertWithCheck(cache, key, result);

    return result;
}


double YB11_lubrication(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_YB::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YB(l)).first;
    }
    const auto& g = g_it->second;

    double ksi = s - 2.0;
    double ksi_inv = 1.0/ksi;
    double log_ksi_inv = std::log(ksi_inv);

    double result = g.g2 * log_ksi_inv + g.g3 * ksi * log_ksi_inv;
    result += B11Y(l, maxIter, rtol, atol);

    return result;
}

double YB11_nf(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_YB::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YB(l)).first;
    }
    const auto& g = g_it->second;
    
    double result = 0.0;
    double two_over_s = 2.0/s;
    double two_over_s_squared = two_over_s * two_over_s;
    double abs_sum = 0.0;
    
    for(int m = 1; m < 2*maxIter; m += 2) {
        double m1 = (m == 2) ? -2.0 : (m-2.0);
        double term = (std::pow(2.0 * (1.0 + l), -m) * fk_YB(m, l) - 
                        2.0/m * g.g2 + 4.0/(m * m1) * g.g3) * std::pow(two_over_s, m);
        
        result += term;
        abs_sum += std::abs(term);
        if(m > 2 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: YB11_nf did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    double log_term = std::log((s+2.0)/(s-2.0));

    result = g.g2 * log_term + g.g3 * (1.0 - two_over_s_squared) * log_term + 2 * g.g3 * two_over_s + result;

    return result;
}


double YB11_ff(double s, double l, int maxIter, double rtol, double atol) {
    double result = 0.0;
    double abs_sum = 0.0;
    double l_plus_1_inv = 1.0/(1.0 + l);
    double s_inv = 1.0/s;
    
    for(int k = 0; k < maxIter; k++) {
        double term = fk_YB(2*k + 1, l) * std::pow(l_plus_1_inv * s_inv, 2 * k + 1);
        result += term;
        abs_sum += std::abs(term);
        
        if(k > 0 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (k == maxIter - 1) {
            std::cerr << "Warning: YB11_ff did not converge within " << maxIter << " iterations." << std::endl;
        }

    }
    
    return result;
}

// ---------------------- YB12 functions ----------------------

/**
 * @brief Computes the B12Y coefficient for the two-body problem.
 * 
 * @details
 * This function calculates the B12Y coefficient using a series expansion with cached values
 * for optimization. It utilizes g-coefficients and implements a convergence check based on
 * relative and absolute tolerances. See Townsend (2023) for the formula.
 *
 * @param l The ratio of the two sphere radii
 * @param maxIter Maximum number of iterations for the series expansion
 * @param rtol Relative tolerance for convergence check
 * @param atol Absolute tolerance for convergence check
 * 
 * @return The computed B12Y coefficient value
 * 
 * @note The function uses caching mechanisms through Cache_YB for both B12Y values
 *       and g-coefficients to improve performance on repeated calls
 * @warning May print a warning to stderr if convergence is not reached within maxIter iterations
 * 
 * @see GCoeffs_YB
 * @see fk_YB
 * @see Cache_YB
 */
double B12Y(double l, int maxIter, double rtol, double atol) {
    // Look for cached value
    auto& cache = Cache_YB::getInstance().BX_cache;
    auto key = std::make_pair(false, l);
    auto it = cache.find(key);
    if(it != cache.end()) {
        return it->second;
    }

    // Find or compute g-coefficients
    auto& g_cache = Cache_YB::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YB(l)).first;
    }
    const auto& g = g_it->second;

    double result = 0.0;
    double abs_sum = 0.0;

    for(int m = 2; m < 2*maxIter; m += 2) {
        double term = (std::pow(2.0 * (1.0 + l), -m) * fk_YB(m, l) - 
                          2.0/m * g.g2 + 4.0/(m * (m + 2.0)) * g.g3);

        result += term;
        abs_sum += std::abs(term);

        if(m > 1 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: B12Y did not converge within " << maxIter << " iterations." << std::endl;
        }

    }

    result = -g.g3 + result;

    // We ommit the prefactor term in the original formula because it's included in the YB12_lubrication function

    // Store result in cache
    Cache_YB::getInstance().insertWithCheck(cache, key, result);

    return result;
}

double YB12_lubrication(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_YB::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YB(l)).first;
    }
    const auto& g = g_it->second;

    double ksi = s - 2.0;
    double ksi_inv = 1.0/ksi;
    double log_ksi_inv = std::log(ksi_inv);

    double result = g.g2 * log_ksi_inv + g.g3 * ksi * log_ksi_inv;
    result += B12Y(l, maxIter, rtol, atol);

    return -result;
}


double YB12_nf(double s, double l, int maxIter, double rtol, double atol) {
    // Get or compute g-coefficients
    auto& g_cache = Cache_YB::getInstance().g_coeffs_cache;
    auto g_it = g_cache.find(l);
    if(g_it == g_cache.end()) {
        g_it = g_cache.emplace(l, GCoeffs_YB(l)).first;
    }
    const auto& g = g_it->second;
    
    double result = 0.0;
    double two_over_s = 2.0/s;
    double two_over_s_squared = two_over_s * two_over_s;
    double abs_sum = 0.0;
    
    for(int m = 2; m < 2*maxIter; m += 2) {
        double m1 = (m == 2) ? -2.0 : (m-2.0);
        double term = (std::pow(2.0 * (1.0 + l), -m) * fk_YB(m, l) - 
                        2.0/m * g.g2 + 4.0/(m * m1) * g.g3) * std::pow(two_over_s, m);
        
        result += term;
        abs_sum += std::abs(term);
        if(m > 1 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (m == 2*maxIter - 2) {
            std::cerr << "Warning: YB12_nf did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    double arg = 1.0 - two_over_s_squared;
    double log_term = std::log((arg));
    result = -g.g2 * log_term - g.g3 * arg * log_term + result;
    
    return -result;
}

double YB12_ff(double s, double l, int maxIter, double rtol, double atol) {
    double result = 0.0;
    double abs_sum = 0.0;
    double l_plus_1_inv = 1.0/(1.0 + l);
    double s_inv = 1.0/s;
    
    for(int k = 0; k < maxIter; k++) {
        double term = fk_YB(2*k, l) * std::pow(l_plus_1_inv * s_inv, 2.0*k);
        result += term;
        abs_sum += std::abs(term);
        
        if(k > 0 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (k == maxIter - 1) {
            std::cerr << "Warning: YB12_ff did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    return -result;
}


// ---------------------- Main functions ----------------------
double YB11(double s, double l, 
            double lubr_cutoff = DEFAULT_LUBR_CUTOFF,
            double cutoff = DEFAULT_CUTOFF,
            int maxIter = DEFAULT_MAX_ITER,
            double rtol = DEFAULT_RTOL,
            double atol = DEFAULT_ATOL) {
    
    if (s < lubr_cutoff) {
        return YB11_lubrication(s, l, maxIter, rtol, atol);
    } else if (s < cutoff) {
        return YB11_nf(s, l, maxIter, rtol, atol);
    } else {
        return YB11_ff(s, l, maxIter, rtol, atol);
    }
}

double YB12(double s, double l, 
            double lubr_cutoff = DEFAULT_LUBR_CUTOFF,
            double cutoff = DEFAULT_CUTOFF,
            int maxIter = DEFAULT_MAX_ITER,
            double rtol = DEFAULT_RTOL,
            double atol = DEFAULT_ATOL) {
    
    if (s < lubr_cutoff) {
        return YB12_lubrication(s, l, maxIter, rtol, atol);
    } else if (s < cutoff) {
        return YB12_nf(s, l, maxIter, rtol, atol);
    } else {
        return YB12_ff(s, l, maxIter, rtol, atol);
    }
}

// Add these implementations before the PYBIND11_MODULE:
namespace yb_utils {
    void clear_cache() {
        Cache_YB::getInstance().clear();
#ifdef ENABLE_PROFILING
        Profiler::getData().clear();
#endif
    }
}
