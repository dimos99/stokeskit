#include "common.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <zlib.h>

struct GlobalCache {
    static const size_t MAX_CACHE_SIZE = 100000000000; // Maximum entries per cache

    std::map<std::pair<int, int>, double> comb_cache;
    std::map<std::pair<double, double>, double> zeta_cache;
    
    static GlobalCache& getInstance() {
        static GlobalCache instance;
        return instance;
    }

#ifdef ENABLE_PROFILING
    void updateMemoryStats() {
        size_t total_bytes = 0;
        size_t total_entries = 0;

        // Calculate memory for comb_cache
        total_bytes += sizeof(std::pair<int, int>) * comb_cache.size() * 2;
        total_bytes += sizeof(double) * comb_cache.size();
        total_entries += comb_cache.size();

        // Calculate memory for zeta_cache
        total_bytes += sizeof(std::pair<double, double>) * zeta_cache.size() * 2;
        total_bytes += sizeof(double) * zeta_cache.size();
        total_entries += zeta_cache.size();

        MemoryProfiler::record("GlobalCache", total_bytes, total_entries);
    }
#endif

    void clear() {
        comb_cache.clear();
        zeta_cache.clear();
#ifdef ENABLE_PROFILING
        updateMemoryStats();
#endif
    }

    template<typename MapType>
    void checkAndTrimCache(MapType& cache) {
        if (cache.size() > MAX_CACHE_SIZE) {
            std::cerr << "Warning: Global cache size exceeded maximum limit. Clearing cache..." << std::endl;
            // Remove last half of elements
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
    GlobalCache() { 
#ifdef ENABLE_PROFILING
        updateMemoryStats(); 
#endif
    } // Private constructor for singleton
};


double comb(int n, int k) {
    auto& cache = GlobalCache::getInstance().comb_cache;
    auto key = std::make_pair(n, k);
    auto it = cache.find(key);
    if(it != cache.end()) {
        return it->second;
    }
    if(k > n || k < 0) return 0.0;
    double r = 1.0;
    for(int i = 1; i <= k; ++i) {
        r = r * (n - k + i) / i;
    }
    GlobalCache::getInstance().insertWithCheck(cache, key, r);
    return r;
}

double zeta_func(double z, double a, int maxIter, double rtol, double atol) {
    auto& cache = GlobalCache::getInstance().zeta_cache;
    auto key = std::make_pair(z, a);
    auto it = cache.find(key);
    if(it != cache.end()) {
        return it->second;
    }

    double result = 0.0;
    double abs_sum = 0.0;

    for(int k = 0; k < maxIter; k++) {
        double term = 1.0 / std::pow(k + a, z);
        result += term;

        abs_sum += std::abs(term);
        if(k > 0 && std::abs(term) < rtol * abs_sum + atol) {
            break;
        }

        // Check if we reached the maximum number of iterations
        if (k == maxIter - 1) {
            std::cerr << "Warning: Zeta function did not converge within " << maxIter << " iterations." << std::endl;
        }

    }

    GlobalCache::getInstance().insertWithCheck(cache, key, result);
    return result;
}


#ifdef ENABLE_PROFILING
// Note: ProfileData is defined in the header as a struct, we don't need to redefine it here

// Profiler method implementations
std::unordered_map<std::string, ProfileData>& Profiler::getData() {
    static std::unordered_map<std::string, ProfileData> data;
    return data;
}

std::mutex& Profiler::getMutex() {
    static std::mutex mutex;
    return mutex;
}

void Profiler::record(const std::string& function, double duration, bool cache_hit) {
    std::lock_guard<std::mutex> lock(getMutex());
    auto& entry = getData()[function];
    entry.calls++;
    entry.total_time += duration;
    if (cache_hit) {
        entry.cache_hits++;
    } else {
        entry.cache_misses++;
    }
}

void Profiler::clear() {
    std::lock_guard<std::mutex> lock(getMutex());
    getData().clear();
}

void Profiler::print_stats() {
    std::lock_guard<std::mutex> lock(getMutex());
    std::cout << "\nPerformance Statistics:\n"
              << "--------------------------------------------------------------------------------\n"
              << std::left << std::setw(20) << "Function"
              << std::setw(10) << "Calls"
              << std::setw(15) << "Total (ms)"
              << std::setw(15) << "Avg (ms)"
              << std::setw(12) << "Cache Hits"
              << "Miss Rate\n"
              << "--------------------------------------------------------------------------------\n";
    for (const auto& entry : getData()) {
        double avg_time = entry.second.calls > 0 
                       ? entry.second.total_time / entry.second.calls 
                       : 0;
        double miss_rate = entry.second.calls > 0 
                        ? static_cast<double>(entry.second.cache_misses) / entry.second.calls 
                        : 0;
        std::cout << std::left << std::setw(20) << entry.first
                  << std::setw(10) << entry.second.calls
                  << std::setw(15) << entry.second.total_time
                  << std::setw(15) << avg_time
                  << std::setw(12) << entry.second.cache_hits
                  << (miss_rate * 100) << "%\n";
    }
}

void MemoryProfiler::print_stats() {
    std::lock_guard<std::mutex> lock(getMutex());
    std::cout << "\nMemory Usage Statistics:\n"
              << "--------------------------------------------------------------------------------\n"
              << std::left << std::setw(30) << "Cache Component"
              << std::setw(15) << "Current MB"
              << std::setw(15) << "Peak MB"
              << std::setw(15) << "Entries"
              << "Peak Entries\n"
              << "--------------------------------------------------------------------------------\n";
    
    double total_mb = 0.0;
    double total_peak_mb = 0.0;
    size_t total_entries = 0;
    size_t total_peak_entries = 0;

    // Group stats by main cache type
    std::map<std::string, std::vector<std::pair<std::string, const MemoryStats*>>> grouped_stats;
    
    for (const auto& entry : getData()) {
        std::string prefix = entry.first.substr(0, entry.first.find('_'));
        std::string component = entry.first;
        grouped_stats[prefix].push_back({component, &entry.second});
    }

    // Print stats by group
    for (const auto& group : grouped_stats) {
        double group_mb = 0.0;
        double group_peak_mb = 0.0;
        size_t group_entries = 0;
        size_t group_peak_entries = 0;

        // Print individual components
        for (const auto& component : group.second) {
            double current_mb = static_cast<double>(component.second->allocated_bytes) / (1024 * 1024);
            double peak_mb = static_cast<double>(component.second->peak_bytes) / (1024 * 1024);
            
            std::cout << std::left << std::setw(30) << component.first
                      << std::setw(15) << std::fixed << std::setprecision(2) << current_mb
                      << std::setw(15) << peak_mb
                      << std::setw(15) << component.second->entry_count
                      << component.second->peak_entries << "\n";

            group_mb += current_mb;
            group_peak_mb += peak_mb;
            group_entries += component.second->entry_count;
            group_peak_entries += component.second->peak_entries;
        }

        // Print group subtotal
        std::cout << "--------------------------------------------------------------------------------\n"
                  << std::left << std::setw(30) << (group.first + " Total")
                  << std::setw(15) << group_mb
                  << std::setw(15) << group_peak_mb
                  << std::setw(15) << group_entries
                  << group_peak_entries << "\n"
                  << "--------------------------------------------------------------------------------\n";

        total_mb += group_mb;
        total_peak_mb += group_peak_mb;
        total_entries += group_entries;
        total_peak_entries += group_peak_entries;
    }

    // Print grand total
    std::cout << std::left << std::setw(30) << "GRAND TOTAL"
              << std::setw(15) << total_mb
              << std::setw(15) << total_peak_mb
              << std::setw(15) << total_entries
              << total_peak_entries << "\n"
              << "--------------------------------------------------------------------------------\n";
}
#endif

namespace cache_utils {

inline size_t npq_to_index(int n, int p, int q, size_t max_p, size_t max_q) {
    return ((size_t)n * max_p * max_q) + ((size_t)p * max_q) + (size_t)q;
}

void ensure_cache_dir(const std::string& dir) {
    if (!std::filesystem::exists(dir)) {
        std::filesystem::create_directory(dir);
    }
}

// Add compression functions
void compress_data(const char* data, size_t size, std::vector<char>& compressed) {
    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    deflateInit(&strm, Z_BEST_COMPRESSION);

    strm.avail_in = size;
    strm.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(data));

    compressed.resize(size);  // Initial size guess
    size_t total_out = 0;

    do {
        if (total_out == compressed.size()) {
            compressed.resize(compressed.size() * 2);
        }
        
        strm.avail_out = compressed.size() - total_out;
        strm.next_out = reinterpret_cast<Bytef*>(compressed.data() + total_out);
        
        deflate(&strm, Z_FINISH);
        total_out = strm.total_out;
    } while (strm.avail_out == 0);

    compressed.resize(total_out);
    deflateEnd(&strm);
}

void decompress_data(const char* compressed_data, size_t compressed_size, std::vector<char>& decompressed) {
    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    inflateInit(&strm);

    strm.avail_in = compressed_size;
    strm.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(compressed_data));

    decompressed.resize(compressed_size * 2);  // Initial size guess
    size_t total_out = 0;

    do {
        if (total_out == decompressed.size()) {
            decompressed.resize(decompressed.size() * 2);
        }
        
        strm.avail_out = decompressed.size() - total_out;
        strm.next_out = reinterpret_cast<Bytef*>(decompressed.data() + total_out);
        
        int ret = inflate(&strm, Z_NO_FLUSH);
        if (ret == Z_STREAM_END) break;
        if (ret != Z_OK) {
            inflateEnd(&strm);
            throw std::runtime_error("Decompression failed");
        }
        total_out = strm.total_out;
    } while (strm.avail_out == 0);

    decompressed.resize(total_out);
    inflateEnd(&strm);
}

void read_cache_file(const std::string& filename, 
                    std::vector<std::pair<bool, double>>& cache,
                    size_t max_n, size_t max_p, size_t max_q) {
    auto start = std::chrono::high_resolution_clock::now();
    
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Warning: Could not open cache file for reading: " << filename << std::endl;
        return;
    }

    std::cout << "Reading cache from: " << filename << std::endl;

    // Read compressed size
    uint32_t compressed_size;
    if (!file.read(reinterpret_cast<char*>(&compressed_size), sizeof(compressed_size))) {
        std::cerr << "Error: Failed to read compressed size from cache file" << std::endl;
        return;
    }

    std::cout << "Reading " << compressed_size << " bytes of compressed data..." << std::endl;

    // Read compressed data
    std::vector<char> compressed_data(compressed_size);
    if (!file.read(compressed_data.data(), compressed_size)) {
        std::cerr << "Error: Failed to read compressed data from cache file" << std::endl;
        return;
    }

    get_stats().compressed_size = compressed_size;
    
    // Decompress the data
    std::vector<char> decompressed_data;
    try {
        std::cout << "Decompressing cache data..." << std::endl;
        decompress_data(compressed_data.data(), compressed_size, decompressed_data);
    } catch (const std::runtime_error& e) {
        std::cerr << "Error during decompression: " << e.what() << std::endl;
        return;
    }

    get_stats().original_size = decompressed_data.size();
    
    // Read header
    const char* current = decompressed_data.data();
    uint32_t magic;
    std::memcpy(&magic, current, sizeof(magic));
    current += sizeof(magic);
    
    if (magic != CACHE_MAGIC) {
        std::cerr << "Error: Invalid cache file format (magic number mismatch) in " << filename << std::endl;
        return;
    }

    // Read entries
    size_t entry_count = 0;
    size_t skipped_entries = 0;
    
    std::cout << "Loading cache entries..." << std::endl;

    while (current < decompressed_data.data() + decompressed_data.size()) {
        uint32_t n, p, q;
        double value;

        std::memcpy(&n, current, sizeof(n));
        current += sizeof(n);
        std::memcpy(&p, current, sizeof(p));
        current += sizeof(p);
        std::memcpy(&q, current, sizeof(q));
        current += sizeof(q);
        std::memcpy(&value, current, sizeof(value));
        current += sizeof(value);

        if (n < max_n && p < max_p && q < max_q) {
            size_t idx = npq_to_index(n, p, q, max_p, max_q);
            cache[idx] = {true, value};
            entry_count++;
        } else {
            skipped_entries++;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    get_stats().read_time = duration;
    get_stats().entries_loaded = entry_count;
    get_stats().entries_skipped = skipped_entries;

    std::cout << "Cache load complete:\n"
              << "  - Loaded entries: " << entry_count << "\n"
              << "  - Skipped entries: " << skipped_entries << "\n"
              << "  - Time taken: " << duration << " ms" << std::endl;
}

void write_cache_file(const std::string& filename,
                     const std::vector<std::pair<bool, double>>& cache,
                     size_t max_n, size_t max_p, size_t max_q) {
    // Skip writing if no modifications were made
    if (!get_stats().cache_modified) {
        std::cout << "Cache unchanged, skipping write operation." << std::endl;
        return;
    }

    auto start = std::chrono::high_resolution_clock::now();
    static bool is_first_write = true;

    std::cout << "Writing cache to: " << filename << std::endl;

    std::ofstream file(filename, std::ios::binary | std::ios::trunc);
    if (!file) {
        std::cerr << "Error: Could not open cache file for writing: " << filename << std::endl;
        return;
    }

    // Count valid entries first
    size_t valid_entries = 0;
    for (const auto& entry : cache) {
        if (entry.first) valid_entries++;
    }

    std::cout << "Preparing to write " << valid_entries << " cache entries..." << std::endl;

    // Prepare the data to be compressed
    std::vector<char> raw_data;
    raw_data.reserve(sizeof(uint32_t) + valid_entries * (4 * sizeof(uint32_t)));

    // Write header
    uint32_t magic = CACHE_MAGIC;
    raw_data.insert(raw_data.end(), 
                   reinterpret_cast<char*>(&magic),
                   reinterpret_cast<char*>(&magic) + sizeof(magic));

    // Write entries
    size_t written_entries = 0;
    for (size_t idx = 0; idx < cache.size(); ++idx) {
        if (cache[idx].first) {
            uint32_t n = idx / (max_p * max_q);
            uint32_t p = (idx / max_q) % max_p;
            uint32_t q = idx % max_q;
            double value = cache[idx].second;

            raw_data.insert(raw_data.end(),
                          reinterpret_cast<const char*>(&n),
                          reinterpret_cast<const char*>(&n) + sizeof(n));
            raw_data.insert(raw_data.end(),
                          reinterpret_cast<const char*>(&p),
                          reinterpret_cast<const char*>(&p) + sizeof(p));
            raw_data.insert(raw_data.end(),
                          reinterpret_cast<const char*>(&q),
                          reinterpret_cast<const char*>(&q) + sizeof(q));
            raw_data.insert(raw_data.end(),
                          reinterpret_cast<const char*>(&value),
                          reinterpret_cast<const char*>(&value) + sizeof(value));
            written_entries++;
        }
    }

    get_stats().original_size = raw_data.size();
    std::cout << "Compressing " << raw_data.size() << " bytes of cache data..." << std::endl;

    // Compress the data
    std::vector<char> compressed_data;
    try {
        compress_data(raw_data.data(), raw_data.size(), compressed_data);
    } catch (const std::exception& e) {
        std::cerr << "Error during compression: " << e.what() << std::endl;
        return;
    }

    get_stats().compressed_size = compressed_data.size();
    get_stats().entries_written = written_entries;

    // Write compressed size and data
    uint32_t compressed_size = compressed_data.size();
    file.write(reinterpret_cast<const char*>(&compressed_size), sizeof(compressed_size));
    if (!file) {
        std::cerr << "Error: Failed to write compressed size" << std::endl;
        return;
    }

    file.write(compressed_data.data(), compressed_data.size());
    if (!file) {
        std::cerr << "Error: Failed to write compressed data" << std::endl;
        return;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();

    if (is_first_write) {
        get_stats().startup_write_time = duration;
        is_first_write = false;
        print_initial_stats();
    } else {
        get_stats().write_time = duration;
    }

    // Reset the modification flag after successful write
    get_stats().cache_modified = false;

    std::cout << "Cache write complete:\n"
              << "  - Written entries: " << written_entries << "\n"
              << "  - Compressed size: " << compressed_size << " bytes\n"
              << "  - Time taken: " << duration << " ms" << std::endl;
}

// Add stats storage

CacheIOStats& get_stats() {
    static CacheIOStats stats; // Will use default member initializers
    return stats;
}

void print_io_stats() {
    const auto& stats = get_stats();
    double compression_ratio = stats.original_size > 0 ? 
        (1.0 - static_cast<double>(stats.compressed_size) / stats.original_size) * 100.0 : 0.0;

    std::cout << "\nCache I/O Statistics:"
              << "\n----------------------------------------"
              << "\nTiming:"
              << "\n  Read time:         " << stats.read_time << " ms"
              << "\n  Write time:        " << stats.write_time << " ms"
              << "\n\nSize Information:"
              << "\n  Original size:     " << stats.original_size << " bytes"
              << "\n  Compressed size:   " << stats.compressed_size << " bytes"
              << "\n  Compression ratio: " << std::fixed << std::setprecision(2) 
              << compression_ratio << "%"
              << "\n\nEntry Statistics:"
              << "\n  Entries loaded:    " << stats.entries_loaded
              << "\n  Entries skipped:   " << stats.entries_skipped
              << "\n  Entries written:   " << stats.entries_written
              << "\n----------------------------------------\n";
}

void print_initial_stats() {
    const auto& stats = get_stats();
    std::cout << "\nCache Initial Write Statistics:\n"
              << "----------------------------------------\n"
              << "Write time:        " << stats.startup_write_time << " ms\n"
              << "----------------------------------------\n";
}

} // namespace cache_utils