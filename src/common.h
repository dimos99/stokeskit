#ifndef COMMON_H
#define COMMON_H

#pragma once

#include <string>
#include <unordered_map>
#include <map>
#include <mutex>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <cmath>    // Add for std::pow
#include <cstring>  // Add for std::memcpy

// Define default parameters as constants
namespace {
    constexpr double DEFAULT_LUBR_CUTOFF = 2.001;
    constexpr double DEFAULT_CUTOFF = 4.0;
    constexpr int DEFAULT_MAX_ITER = 200;
    constexpr double DEFAULT_RTOL = 1e-4;
    constexpr double DEFAULT_ATOL = 1e-6;
}

// Cache-related declarations
namespace cache_utils {
    constexpr uint32_t CACHE_MAGIC = 0x43414348; // 'CACH' in hex
    constexpr size_t CHUNK = 16384;  // Size of compression buffer chunks
    
    void ensure_cache_dir(const std::string& dir);
    
    void read_cache_file(const std::string& filename, 
                        std::vector<std::pair<bool, double>>& cache,
                        size_t max_n, size_t max_p, size_t max_q);
    
    void write_cache_file(const std::string& filename,
                         const std::vector<std::pair<bool, double>>& cache,
                         size_t max_n, size_t max_p, size_t max_q);
    
    void compress_data(const char* data, size_t size, std::vector<char>& compressed);
    void decompress_data(const char* compressed_data, size_t compressed_size, std::vector<char>& decompressed);

    struct CacheIOStats {
        double read_time = 0.0;
        double write_time = 0.0;
        double startup_write_time = 0.0;
        size_t original_size = 0;
        size_t compressed_size = 0;
        size_t entries_loaded = 0;
        size_t entries_skipped = 0;
        size_t entries_written = 0;
        bool cache_modified = false;  // Add this line
    };

    CacheIOStats& get_stats();
    void print_io_stats();
    void print_initial_stats();  // Add this new function declaration
}

// Helper function for binomial coefficients (nCk).
/**
 * @brief Calculates the binomial coefficient C(n,k) using dynamic programming
 * 
 * This function computes the combination of n things taken k at a time using
 * a cache-based approach for improved performance. The result is stored in
 * a global cache to avoid redundant calculations.
 * 
 * @param n The total number of elements
 * @param k The number of elements to choose
 * @return double The binomial coefficient C(n,k)
 * 
 * @note The function uses a cache to store previously calculated values
 * @note Returns 0 if k > n or k < 0
 * 
 * Time Complexity: O(k) for uncached values, O(1) for cached values
 * Space Complexity: O(1) for the calculation, O(n*k) for the cache
 */
double comb(int n, int k);

// Helper function for zeta function (Hurwitz zeta function).
/**
 * @brief Calculates the Hurwitz zeta function for a given complex number
 * 
 * This function computes the Hurwitz zeta function for a given complex number
 * using a series expansion approach. The function uses a combination of
 * logarithmic terms and zeta function evaluations to compute the result.
 * 
 * @param s The complex number parameter
 * @param a The real part of the complex number
 * @param maxIter Maximum number of iterations for convergence
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * @return double The value of the Hurwitz zeta function
 * 
 * @note The function uses a series expansion for the zeta function
 * @note The function uses a combination of logarithmic terms and zeta evaluations
 * 
 * Time Complexity: O(maxIter) for convergence
 * Space Complexity: O(1)
 */
double zeta_func(double z, double a, int maxIter, double rtol, double atol);


// Common profiling functionality
#ifdef ENABLE_PROFILING
struct ProfileData {
    size_t calls = 0;
    double total_time = 0.0;
    size_t cache_hits = 0;
    size_t cache_misses = 0;
};

struct MemoryStats {
    size_t allocated_bytes{0};
    size_t peak_bytes{0};
    size_t entry_count{0};
    size_t peak_entries{0};

    void update(size_t current_bytes, size_t current_entries) {
        allocated_bytes = current_bytes;
        entry_count = current_entries;
        peak_bytes = std::max(peak_bytes, current_bytes);
        peak_entries = std::max(peak_entries, current_entries);
    }
};

class MemoryProfiler {
public:
    static std::unordered_map<std::string, MemoryStats>& getData() {
        static std::unordered_map<std::string, MemoryStats> data;
        return data;
    }

    static std::mutex& getMutex() {
        static std::mutex mutex;
        return mutex;
    }

    static void record(const std::string& cache_name, size_t bytes, size_t entries) {
        std::lock_guard<std::mutex> lock(getMutex());
        getData()[cache_name].update(bytes, entries);
    }

    static void clear() {
        std::lock_guard<std::mutex> lock(getMutex());
        getData().clear();
    }

    static void print_stats();
};

class Profiler {
public:
    static std::unordered_map<std::string, ProfileData>& getData();
    static std::mutex& getMutex();
    static void record(const std::string& function, double duration, bool cache_hit);
    static void print_stats();
    static void clear();
};

// Helper macro for timing blocks of code
#define PROFILE_BLOCK(name) \
    auto start = std::chrono::high_resolution_clock::now(); \
    bool cache_hit = false; \
    auto profile_end = std::make_unique<std::function<void()>>([&](){ \
        auto end = std::chrono::high_resolution_clock::now(); \
        auto duration = std::chrono::duration<double, std::milli>(end - start).count(); \
        Profiler::record(name, duration, cache_hit); \
    }); \

#define PROFILE_CACHE_HIT cache_hit = true
#else
#define PROFILE_BLOCK(name)
#define PROFILE_CACHE_HIT
#endif

#endif // COMMON_H