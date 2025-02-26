#include "XC.h"

#include "common.h"


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
 * - AX_cache: Combined cache for B11Y and B12Y functions, indexed by type flag and l parameter
 * - comb_cache: Stores combinatorial values indexed by integer pairs
 * 
 * The singleton pattern ensures that only one instance of the cache exists throughout
 * the program's lifetime, providing a centralized storage for computed values.
 * 
 * Usage:
 * @code
 * Cache_XC& cache = Cache_XC::getInstance();
 * cache.clear(); // To clear all cached values
 * @endcode
 */
struct Cache_XC {
    static const size_t MAX_CACHE_SIZE = 100000000000; // Updated size to match XA

    std::map<std::tuple<int,int,int>, double> Q_cache; // (n, p, q) -> value
    std::map<std::pair<int, double>, double> fk_cache; // (k, l) -> value

    
    static Cache_XC& getInstance() {
        static Cache_XC instance;
        return instance;
    }
#ifdef ENABLE_PROFILING
    void updateMemoryStats() {
        size_t total_bytes = 0;
        size_t total_entries = 0;

        // Q_cache memory
        total_bytes += sizeof(std::tuple<int,int,int>) * Q_cache.size() * 2;
        total_bytes += sizeof(double) * Q_cache.size();
        total_entries += Q_cache.size();

        // fk_cache memory
        total_bytes += sizeof(std::pair<int, double>) * fk_cache.size() * 2;
        total_bytes += sizeof(double) * fk_cache.size();
        total_entries += fk_cache.size();

        MemoryProfiler::record("Cache_XC", total_bytes, total_entries);
    }
#endif

    void clear() {
        Q_cache.clear();
        fk_cache.clear();
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
    Cache_XC() { 
#ifdef ENABLE_PROFILING
        updateMemoryStats(); 
#endif
    }
};


/** 
 * @brief Calculates the Q_{n,p,q} pure function for the XC scalar function in a two-body resistance problem.
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
 * @note The function uses Cache_XC for memoization
 * 
 */
double Q_npq_XC(int n, int p, int q) {
#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    auto& cache = Cache_XC::getInstance().Q_cache;
    auto key = std::make_tuple(n, p, q);
    auto it = cache.find(key);
    if(it != cache.end()) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("Q_npq_XC", duration, cache_hit);
#endif
        return it->second;
    }

    // Base case
    if(p == 0 && q == 0) {
        double Q = (n == 1) ? 1.0 : 0.0;
        Cache_XC::getInstance().insertWithCheck(cache, key, Q);
#ifdef ENABLE_PROFILING
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("Q_npq_XC", duration, cache_hit);
#endif
        return Q;
    }


    // Recursive formula
    double pre_term = 1.0 / (n + 1.0);
    double Q = 0.0;
    for(int s = 1; s <= q; ++s) {
        Q += comb(n + s, n) * (
            pre_term * s * Q_npq_XC(s, q - s - 1, p - n)
        );
    }

    Cache_XC::getInstance().insertWithCheck(cache, key, Q);

#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("Q_npq_XC", duration, cache_hit);
#endif
    return Q;
}


/**
 * @brief Calculates the fk function for the XC scalar function in a two-body resistance problem.
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
 * @see Q_npq_XC, Cache_XC
 */
double fk_XC(int k, double l) {
#ifdef ENABLE_PROFILING
    auto start = std::chrono::high_resolution_clock::now();
    bool cache_hit = false;
#endif

    auto& cache = Cache_XC::getInstance().fk_cache;
    auto key = std::make_pair(k, l);
    auto it = cache.find(key);
    if(it != cache.end()) {
#ifdef ENABLE_PROFILING
        cache_hit = true;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start).count();
        Profiler::record("fk_XC", duration, cache_hit);
#endif
        return it->second;
    }

    double result = 0.0;
    

    for(int q = 0; q <= k; ++q) {
        // #ifdef _OPENMP
        // std::cout << "[fk_XC debug] Thread " << omp_get_thread_num() 
        //           << " handling q = " << q 
        //           << std::endl;
        // #endif
        result += Q_npq_XC(1, k - q, q) * std::pow(l, q + k%2);
    }
    result *= std::pow(2.0, k);

    Cache_XC::getInstance().insertWithCheck(cache, key, result);
#ifdef ENABLE_PROFILING
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    Profiler::record("fk_XC", duration, cache_hit);
#endif
    return result;
}

// ---------------------- XC11 functions ----------------------

// Lubrication expression

double XC11_lubrication(double s, double l, int maxIter, double rtol, double atol) {

    double ksi = s - 2.0;
    double ksi_inv = 1.0/ksi;
    double log_ksi_inv = std::log(ksi_inv);
    double l2 = l * l;
    double l3 = l2 * l;
    double one_plus_l = 1.0 + l;

    double term1 = l3 / std::pow(one_plus_l, 3);
    double term2 = - l2 / (4 * one_plus_l);


    double result = term1 * zeta_func(3.0, l / one_plus_l, maxIter, rtol, atol) + term2 * ksi * log_ksi_inv;

    return result;
}

double XC11_nf(double s, double l, int maxIter, double rtol, double atol) {

    double result = 0.0;
    double l2 = l * l;
    double one_plus_l = 1.0 + l;
    double one_plus_l_inv = 1.0 / one_plus_l;
    double one_over_s = 1.0 / s;
    double one_over_s_squared = one_over_s * one_over_s;

    double abs_sum = 0.0;
    
    for(int k = 1; k < maxIter; k++) {
        double term = (std::pow(one_plus_l_inv, 2*k) * fk_XC(2*k, l) - 
                        std::pow(2.0, 2*k + 1) / (k * (2*k - 1)) * l2 / 4 * one_plus_l_inv) * 
                        std::pow(one_over_s_squared, k);
        
        result += term;

        abs_sum += std::abs(term);
        if(k > 1 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (k == maxIter - 1) {
            std::cerr << "Warning: XC11_nf did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    double log_term = std::log((s+2.0)/(s-2.0));

    result = l2 / 2 * one_plus_l_inv * std::log(1 - 4 * one_over_s_squared) + 
                l2 * one_plus_l_inv * one_over_s * log_term +
                1 +
                result;

    return result;
}


double XC11_ff(double s, double l, int maxIter, double rtol, double atol) {
    double result = 0.0;
    double abs_sum = 0.0;
    double l_plus_1_inv = 1.0/(1.0 + l);
    double s_inv = 1.0/s;
    
    for(int k = 0; k < maxIter; k++) {
        double term = fk_XC(2*k, l) * std::pow(l_plus_1_inv * s_inv, 2*k);
        result += term;
        abs_sum += std::abs(term);
        
        if(k > 0 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (k == maxIter - 1) {
            std::cerr << "Warning: XC11_ff did not converge within " << maxIter << " iterations." << std::endl;
        }

    }
    
    return result;
}

// ---------------------- XC12 functions ----------------------

// Lubrication expression
double XC12_lubrication(double s, double l, int maxIter, double rtol, double atol) {

    double l2 = l * l;
    double l3 = l2 * l;
    double one_plus_l = 1.0 + l;
    double ksi = s - 2.0;
    double ksi_inv = 1.0/ksi;
    double log_ksi_inv = std::log(ksi_inv);

    double term1 = - l3 / std::pow(one_plus_l, 3);
    double term2 = l2 / (4 * one_plus_l);

    double result = term1 * zeta_func(3.0, 1.0, maxIter, rtol, atol) + term2 * ksi * log_ksi_inv;

    return result;
}


double XC12_nf(double s, double l, int maxIter, double rtol, double atol) {
    double l2 = l * l;
    double one_plus_l = 1.0 + l;
    double one_plus_l_inv = 1.0 / one_plus_l;
    double s_inv = 1.0/s;
    double two_over_s = 2.0/s;
    double two_over_s_squared = two_over_s * two_over_s;
    
    double result = 0.0;
    double abs_sum = 0.0;

    for(int k = 1; k < maxIter; k++) {
        double term = (std::pow(one_plus_l_inv, 2*k + 1) * fk_XC(2*k + 1, l) - 
                        std::pow(2.0, 2*k + 2) / (k * (2*k + 1)) * l2 / 4 * one_plus_l_inv) * 
                        std::pow(s_inv, 2*k + 1);
        
        result += term;
        abs_sum += std::abs(term);
        if(k > 1 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        if (k == maxIter - 1) {
            std::cerr << "Warning: XC12_nf did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    double arg = 1.0 - two_over_s_squared;
    double log_term = std::log((s+2.0)/(s-2.0));
    double l2_over_oneplusl = l2 / one_plus_l;
    
    result = (log_term / 2 + s_inv * std::log(arg) - 2 * s_inv) * l2_over_oneplusl -
                result;
    
    return result;
}

double XC12_ff(double s, double l, int maxIter, double rtol, double atol) {
    double result = 0.0;
    double abs_sum = 0.0;
    double l_plus_1_inv = 1.0/(1.0 + l);
    double s_inv = 1.0/s;
    
    for(int k = 0; k < maxIter; k++) {
        double term = fk_XC(2*k + 1, l) * std::pow(l_plus_1_inv * s_inv, 2.0*k + 1);
        result += term;
        abs_sum += std::abs(term);
        
        if(k > 0 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (k == maxIter - 1) {
            std::cerr << "Warning: XC12_ff did not converge within " << maxIter << " iterations." << std::endl;
        }
    }
    
    return -result;
}


// ---------------------- Main functions ----------------------
double XC11(double s, double l, 
            double lubr_cutoff = DEFAULT_LUBR_CUTOFF,
            double cutoff = DEFAULT_CUTOFF,
            int maxIter = DEFAULT_MAX_ITER,
            double rtol = DEFAULT_RTOL,
            double atol = DEFAULT_ATOL) {
    
    if (s < lubr_cutoff) {
        return XC11_lubrication(s, l, maxIter, rtol, atol);
    } else if (s < cutoff) {
        return XC11_nf(s, l, maxIter, rtol, atol);
    } else {
        return XC11_ff(s, l, maxIter, rtol, atol);
    }
}

double XC12(double s, double l, 
            double lubr_cutoff = DEFAULT_LUBR_CUTOFF,
            double cutoff = DEFAULT_CUTOFF,
            int maxIter = DEFAULT_MAX_ITER,
            double rtol = DEFAULT_RTOL,
            double atol = DEFAULT_ATOL) {
    
    if (s < lubr_cutoff) {
        return XC12_lubrication(s, l, maxIter, rtol, atol);
    } else if (s < cutoff) {
        return XC12_nf(s, l, maxIter, rtol, atol);
    } else {
        return XC12_ff(s, l, maxIter, rtol, atol);
    }
}

// Add these implementations before the PXCIND11_MODULE:
namespace xc_utils {
    void clear_cache() {
        Cache_XC::getInstance().clear();
#ifdef ENABLE_PROFILING
        Profiler::getData().clear();
#endif
    }
}
