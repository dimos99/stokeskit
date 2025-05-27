#ifndef CACHE_HASH_H
#define CACHE_HASH_H

#include <functional>
#include <tuple>
#include <utility>

// Hash functions for complex key types used in caches

// Hash function for std::tuple<int, int, int>
struct tuple_int3_hash {
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
        auto h1 = std::hash<int>{}(std::get<0>(t));
        auto h2 = std::hash<int>{}(std::get<1>(t));
        auto h3 = std::hash<int>{}(std::get<2>(t));
        // Combine hashes using a better mixing function
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

// Hash function for std::pair<int, double>
struct pair_int_double_hash {
    std::size_t operator()(const std::pair<int, double>& p) const {
        auto h1 = std::hash<int>{}(p.first);
        auto h2 = std::hash<double>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};

// Hash function for std::pair<bool, double>
struct pair_bool_double_hash {
    std::size_t operator()(const std::pair<bool, double>& p) const {
        auto h1 = std::hash<bool>{}(p.first);
        auto h2 = std::hash<double>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};

// Hash function for std::pair<int, int>
struct pair_int_int_hash {
    std::size_t operator()(const std::pair<int, int>& p) const {
        auto h1 = std::hash<int>{}(p.first);
        auto h2 = std::hash<int>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};

// Hash function for std::pair<double, double>
struct pair_double_double_hash {
    std::size_t operator()(const std::pair<double, double>& p) const {
        auto h1 = std::hash<double>{}(p.first);
        auto h2 = std::hash<double>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};

#endif // CACHE_HASH_H
