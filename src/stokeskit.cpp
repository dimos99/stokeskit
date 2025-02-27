#include "stokeskit.h"
#include <pybind11/numpy.h>
#include <omp.h>

void init_XA(py::module& m) {
    auto xa = m.def_submodule("XA", "XA scalar resistance functions");
    
    xa.def("P_npq", &P_npq_XA,
        "Calculate P coefficient for XA scalar functions",
        py::arg("n"), py::arg("p"), py::arg("q"));

    xa.def("V_npq", &V_npq_XA,
        "Calculate V coefficient for XA scalar functions",
        py::arg("n"), py::arg("p"), py::arg("q"));

    xa.def("fk", &fk_XA,
        "Calculate fk coefficient for XA scalar functions",
        py::arg("k"), py::arg("l"));

    xa.def("XA11", &XA11,
        "Calculate XA11 scalar resistance function",
        py::arg("s"), py::arg("l"), 
        py::arg("lubr_cutoff")=2.001, // TODO: Decide on default values
        py::arg("cutoff")=4.0,
        py::arg("maxIter")=200,
        py::arg("rtol")=1e-4,
        py::arg("atol")=1e-6);

    xa.def("XA12", &XA12,
        "Calculate XA12 scalar resistance function",
        py::arg("s"), py::arg("l"),
        py::arg("lubr_cutoff")=2.001, // TODO: Decide on default values
        py::arg("cutoff")=4.0,
        py::arg("maxIter")=200,
        py::arg("rtol")=1e-4,
        py::arg("atol")=1e-6);

    xa.def("clear_cache", &xa_utils::clear_cache,
        "Clear all cached values and reset profiling data for XA functions");
}

void init_YA(py::module& m) {
    auto ya = m.def_submodule("YA", "YA scalar resistance functions");
    
    ya.def("P_npq", &P_npq_YA,
        "Calculate P coefficient for YA scalar functions",
        py::arg("n"), py::arg("p"), py::arg("q"));

    ya.def("V_npq", &V_npq_YA,
        "Calculate V coefficient for YA scalar functions",
        py::arg("n"), py::arg("p"), py::arg("q"));

    ya.def("Q_npq", &Q_npq_YA,
        "Calculate Q coefficient for YA scalar functions",
        py::arg("n"), py::arg("p"), py::arg("q"));

    ya.def("fk", &fk_YA,
        "Calculate fk coefficient for YA scalar functions",
        py::arg("k"), py::arg("l"));

    ya.def("YA11", &YA11,
        "Calculate XA11 scalar resistance function",
        py::arg("s"), py::arg("l"), 
        py::arg("lubr_cutoff")=2.001, // TODO: Decide on default values
        py::arg("cutoff")=4.0,
        py::arg("maxIter")=200,
        py::arg("rtol")=1e-4,
        py::arg("atol")=1e-6);

    ya.def("YA12", &YA12,
        "Calculate XA12 scalar resistance function",
        py::arg("s"), py::arg("l"),
        py::arg("lubr_cutoff")=2.001, // TODO: Decide on default values
        py::arg("cutoff")=4.0,
        py::arg("maxIter")=200,
        py::arg("rtol")=1e-4,
        py::arg("atol")=1e-6);

    ya.def("clear_cache", &ya_utils::clear_cache,
        "Clear all cached values and reset profiling data for XA functions");
}

void init_YB(py::module& m) {
    auto yb = m.def_submodule("YB", "YB scalar resistance functions");
    
    yb.def("P_npq", &P_npq_YB,
        "Calculate P coefficient for YB scalar functions",
        py::arg("n"), py::arg("p"), py::arg("q"));

    yb.def("V_npq", &V_npq_YB,
        "Calculate V coefficient for YB scalar functions",
        py::arg("n"), py::arg("p"), py::arg("q"));

    yb.def("Q_npq", &Q_npq_YB,
        "Calculate Q coefficient for YB scalar functions",
        py::arg("n"), py::arg("p"), py::arg("q"));

    yb.def("fk", &fk_YB,
        "Calculate fk coefficient for YB scalar functions",
        py::arg("k"), py::arg("l"));

    yb.def("YB11", &YB11,
        "Calculate XA11 scalar resistance function",
        py::arg("s"), py::arg("l"), 
        py::arg("lubr_cutoff")=2.001, // TODO: Decide on default values
        py::arg("cutoff")=4.0,
        py::arg("maxIter")=200,
        py::arg("rtol")=1e-4,
        py::arg("atol")=1e-6);

    yb.def("YB12", &YB12,
        "Calculate XA12 scalar resistance function",
        py::arg("s"), py::arg("l"),
        py::arg("lubr_cutoff")=2.001, // TODO: Decide on default values
        py::arg("cutoff")=4.0,
        py::arg("maxIter")=200,
        py::arg("rtol")=1e-4,
        py::arg("atol")=1e-6);

    yb.def("clear_cache", &ya_utils::clear_cache,
        "Clear all cached values and reset profiling data for XA functions");
}


void init_XC(py::module& m) {
    auto xc = m.def_submodule("XC", "XC scalar resistance functions");
    
    xc.def("Q_npq", &Q_npq_XC,
        "Calculate Q coefficient for XC scalar functions",
        py::arg("n"), py::arg("p"), py::arg("q"));

    xc.def("fk", &fk_XC,
        "Calculate fk coefficient for XC scalar functions",
        py::arg("k"), py::arg("l"));

    xc.def("XC11", &XC11,
        "Calculate XC11 scalar resistance function",
        py::arg("s"), py::arg("l"), 
        py::arg("lubr_cutoff")=2.001, // TODO: Decide on default values
        py::arg("cutoff")=4.0,
        py::arg("maxIter")=200,
        py::arg("rtol")=1e-4,
        py::arg("atol")=1e-6);

    xc.def("XC12", &XC12,
        "Calculate XC12 scalar resistance function",
        py::arg("s"), py::arg("l"),
        py::arg("lubr_cutoff")=2.001, // TODO: Decide on default values
        py::arg("cutoff")=4.0,
        py::arg("maxIter")=200,
        py::arg("rtol")=1e-4,
        py::arg("atol")=1e-6);

    xc.def("clear_cache", &xc_utils::clear_cache,
        "Clear all cached values and reset profiling data for XC functions");
}

void init_YC(py::module& m) {
    auto yc = m.def_submodule("YC", "YC scalar resistance functions");
    
    yc.def("Q_npq", &Q_npq_YC,
        "Calculate Q coefficient for YC scalar functions",
        py::arg("n"), py::arg("p"), py::arg("q"));

    yc.def("fk", &fk_YC,
        "Calculate fk coefficient for YC scalar functions",
        py::arg("k"), py::arg("l"));

    yc.def("YC11", &YC11,
        "Calculate YC11 scalar resistance function",
        py::arg("s"), py::arg("l"), 
        py::arg("lubr_cutoff")=2.001, // TODO: Decide on default values
        py::arg("cutoff")=4.0,
        py::arg("maxIter")=200,
        py::arg("rtol")=1e-4,
        py::arg("atol")=1e-6);

    yc.def("YC12", &YC12,
        "Calculate YC12 scalar resistance function",
        py::arg("s"), py::arg("l"),
        py::arg("lubr_cutoff")=2.001, // TODO: Decide on default values
        py::arg("cutoff")=4.0,
        py::arg("maxIter")=200,
        py::arg("rtol")=1e-4,
        py::arg("atol")=1e-6);

    yc.def("clear_cache", &yc_utils::clear_cache,
        "Clear all cached values and reset profiling data for YC functions");
}

void init_Minfinity(py::module& m) {
    auto minf = m.def_submodule("Minfinity", "Mobility matrices for many-body hydrodynamic interactions");
    
    // Expose Levi-Civita tensor
    minf.def("levi", &levi, "Levi-Civita tensor");

    // Expose tensor operations
    minf.def("J_tensor", &J_tensor, "J tensor calculation");
    minf.def("Lap_J", &Lap_J, "Laplacian of J tensor");
    minf.def("D_J", &D_J, "Derivative of J tensor");
    minf.def("R_tensor", &R_tensor, "R tensor calculation");
    minf.def("K_tensor", &K_tensor, "K tensor calculation");
    minf.def("DD_J", &DD_J, "Second derivative of J tensor");
    minf.def("DLap_J", &DLap_J, "Derivative of Laplacian of J tensor");
    
    // Expose main computation functions - updated implementation
    minf.def("computeMinfinityArray", [](const py::array_t<double>& positions_array, 
                                      const py::array_t<double>& radii_array, double mu) {
        // Check input arrays
        if (positions_array.ndim() != 2 || positions_array.shape(1) != 3) {
            throw std::runtime_error("positions must be an Nx3 array");
        }
        if (radii_array.ndim() != 1) {
            throw std::runtime_error("radii must be a 1D array");
        }
        if (positions_array.shape(0) != radii_array.shape(0)) {
            throw std::runtime_error("Number of positions must match number of radii");
        }

        // Convert numpy arrays to vectors
        std::vector<Eigen::Vector3d> positions;
        std::vector<double> radii;
        
        for (py::ssize_t i = 0; i < positions_array.shape(0); i++) {
            Eigen::Vector3d pos(
                positions_array.at(i, 0),
                positions_array.at(i, 1),
                positions_array.at(i, 2)
            );
            positions.push_back(pos);
            radii.push_back(radii_array.at(i));
        }
        
        // Call the C++ function
        return computeMinfinity(positions, radii, mu);
    }, "Compute the Minfinity mobility matrix for a collection of particles");
    
    // Add Rinfinity computation
    minf.def("computeRinfinityArray", [](const py::array_t<double>& positions_array, 
                                      const py::array_t<double>& radii_array, double mu) {
        // Check input arrays
        if (positions_array.ndim() != 2 || positions_array.shape(1) != 3) {
            throw std::runtime_error("positions must be an Nx3 array");
        }
        if (radii_array.ndim() != 1) {
            throw std::runtime_error("radii must be a 1D array");
        }
        if (positions_array.shape(0) != radii_array.shape(0)) {
            throw std::runtime_error("Number of positions must match number of radii");
        }

        // Convert numpy arrays to vectors
        std::vector<Eigen::Vector3d> positions;
        std::vector<double> radii;
        
        for (py::ssize_t i = 0; i < positions_array.shape(0); i++) {
            Eigen::Vector3d pos(
                positions_array.at(i, 0),
                positions_array.at(i, 1),
                positions_array.at(i, 2)
            );
            positions.push_back(pos);
            radii.push_back(radii_array.at(i));
        }
        
        // Call the C++ function
        return computeRinfinity(positions, radii, mu);
    }, "Compute the Rinfinity resistance matrix (inverse of Minfinity) for a collection of particles");

    // Expose OpenMP-related functions if available
    #ifdef _OPENMP
    minf.def("get_max_threads", []() { 
        return omp_get_max_threads(); 
    }, "Get maximum number of OpenMP threads available");
    
    minf.def("set_num_threads", [](int num_threads) { 
        if (num_threads > 0) {
            omp_set_num_threads(num_threads);
        }
    }, "Set number of OpenMP threads");
    
    minf.def("get_num_threads", []() { 
        int num_threads = 0;
        #pragma omp parallel
        {
            #pragma omp single
            num_threads = omp_get_num_threads();
        }
        return num_threads;
    }, "Get current number of OpenMP threads");
    
    minf.def("has_openmp", []() { return true; }, "Check if OpenMP is available");
    #else
    minf.def("get_max_threads", []() { return 1; }, 
        "Get maximum number of threads available (always 1 without OpenMP)");
    
    minf.def("set_num_threads", [](int num_threads) { /* Do nothing */ }, 
        "Set number of threads (no effect without OpenMP)");
    
    minf.def("get_num_threads", []() { return 1; }, 
        "Get current number of threads (always 1 without OpenMP)");
    
    minf.def("has_openmp", []() { return false; }, "Check if OpenMP is available");
    #endif
}

/**
 * @brief Compute all resistance functions for an array of distances
 * 
 * @param distances Numpy array of distances
 * @param l_values Numpy array of radius ratios
 * @param lubr_cutoff Lubrication cutoff distance
 * @param cutoff Far-field cutoff distance
 * @param maxIter Maximum number of iterations
 * @param rtol Relative tolerance
 * @param atol Absolute tolerance
 * @return py::dict Dictionary containing numpy arrays of computed functions
 */
py::dict compute_resistance_functions(
    py::array_t<double, py::array::c_style | py::array::forcecast> distances,
    py::array_t<double, py::array::c_style | py::array::forcecast> l_values,
    double lubr_cutoff = DEFAULT_LUBR_CUTOFF,
    double cutoff = DEFAULT_CUTOFF,
    int maxIter = DEFAULT_MAX_ITER,
    double rtol = DEFAULT_RTOL,
    double atol = DEFAULT_ATOL) {
    
    if (!distances.size() || !l_values.size()) {
        throw std::runtime_error("Input arrays cannot be empty");
    }

    if (distances.size() != l_values.size()) {
        throw std::runtime_error("distances and l_values arrays must have the same size");
    }

    // Create result arrays filled with zeros
    size_t n = distances.size();
    std::vector<py::array_t<double>> results_arrays;
    results_arrays.reserve(10);  // We know we'll have exactly 10 arrays

    // Pre-allocate all arrays with zeros
    for (int i = 0; i < 10; ++i) {
        results_arrays.push_back(py::array_t<double>(py::ssize_t(n)));
        std::memset(results_arrays.back().mutable_data(), 0, sizeof(double) * n);
    }

    // Get raw pointers for input arrays - these are guaranteed to be contiguous due to forcecast
    const double* s = distances.data();
    const double* l = l_values.data();

    // Get raw pointers for output arrays
    double* xa11_ptr = static_cast<double*>(results_arrays[0].mutable_data());
    double* xa12_ptr = static_cast<double*>(results_arrays[1].mutable_data());
    double* ya11_ptr = static_cast<double*>(results_arrays[2].mutable_data());
    double* ya12_ptr = static_cast<double*>(results_arrays[3].mutable_data());
    double* yb11_ptr = static_cast<double*>(results_arrays[4].mutable_data());
    double* yb12_ptr = static_cast<double*>(results_arrays[5].mutable_data());
    double* xc11_ptr = static_cast<double*>(results_arrays[6].mutable_data());
    double* xc12_ptr = static_cast<double*>(results_arrays[7].mutable_data());
    double* yc11_ptr = static_cast<double*>(results_arrays[8].mutable_data());
    double* yc12_ptr = static_cast<double*>(results_arrays[9].mutable_data());

    // Larger chunk size for better performance
    const size_t CHUNK_SIZE = 5000;
    
    try {
        #pragma omp parallel for schedule(dynamic) 
        for(size_t chunk_start = 0; chunk_start < n; chunk_start += CHUNK_SIZE) {
            size_t chunk_end = std::min(chunk_start + CHUNK_SIZE, n);
            
            for(size_t i = chunk_start; i < chunk_end; i++) {
                // Skip computation if s < 2 or l < 0
                if (s[i] < 2.0 || l[i] < 0.0) {
                    xa11_ptr[i] = 0.0;
                    xa12_ptr[i] = 0.0;
                    ya11_ptr[i] = 0.0;
                    ya12_ptr[i] = 0.0;
                    yb11_ptr[i] = 0.0;
                    yb12_ptr[i] = 0.0;
                    xc11_ptr[i] = 0.0;
                    xc12_ptr[i] = 0.0;
                    yc11_ptr[i] = 0.0;
                    yc12_ptr[i] = 0.0;
                    continue;
                }
                
                xa11_ptr[i] = XA11(s[i], l[i], lubr_cutoff, cutoff, maxIter, rtol, atol);
                xa12_ptr[i] = XA12(s[i], l[i], lubr_cutoff, cutoff, maxIter, rtol, atol);
                ya11_ptr[i] = YA11(s[i], l[i], lubr_cutoff, cutoff, maxIter, rtol, atol);
                ya12_ptr[i] = YA12(s[i], l[i], lubr_cutoff, cutoff, maxIter, rtol, atol);
                yb11_ptr[i] = YB11(s[i], l[i], lubr_cutoff, cutoff, maxIter, rtol, atol);
                yb12_ptr[i] = YB12(s[i], l[i], lubr_cutoff, cutoff, maxIter, rtol, atol);
                xc11_ptr[i] = XC11(s[i], l[i], lubr_cutoff, cutoff, maxIter, rtol, atol);
                xc12_ptr[i] = XC12(s[i], l[i], lubr_cutoff, cutoff, maxIter, rtol, atol);
                yc11_ptr[i] = YC11(s[i], l[i], lubr_cutoff, cutoff, maxIter, rtol, atol);
                yc12_ptr[i] = YC12(s[i], l[i], lubr_cutoff, cutoff, maxIter, rtol, atol);
            }
        }
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Computation failed: ") + e.what());
    }

    // Create results dictionary without moving arrays
    py::dict results;
    results["XA11"] = results_arrays[0];
    results["XA12"] = results_arrays[1];
    results["YA11"] = results_arrays[2];
    results["YA12"] = results_arrays[3];
    results["YB11"] = results_arrays[4];
    results["YB12"] = results_arrays[5];
    results["XC11"] = results_arrays[6];
    results["XC12"] = results_arrays[7];
    results["YC11"] = results_arrays[8];
    results["YC12"] = results_arrays[9];

    return results;
}

PYBIND11_MODULE(stokeskit, m) {
    m.doc() = "Resistance scalar functions for two-body hydrodynamic interactions";
    
#ifdef ENABLE_PROFILING
    // Add memory profiling function
    m.def("print_memory_usage", []() {
        Profiler::print_stats();        // Print performance stats
        MemoryProfiler::print_stats();  // Print memory stats
    }, "Print detailed memory usage statistics for all caches");
#endif

    // Add print profiling stats function
#ifdef ENABLE_PROFILING
    m.def("print_profile_stats", &Profiler::print_stats,
        "Print profiling statistics for all functions");
#endif

    // Add zeta function
    m.def("zeta_func", &zeta_func,
        "Calculate the Hurwitz zeta function for a given complex number",
        py::arg("z"), py::arg("a"), py::arg("maxIter")=200, py::arg("rtol")=1e-4, py::arg("atol")=1e-6);
        
    m.def("compute_resistance_functions", &compute_resistance_functions,
        "Compute all resistance functions for arrays of distances and radius ratios",
        py::arg("distances"),
        py::arg("l_values"),
        py::arg("lubr_cutoff")=DEFAULT_LUBR_CUTOFF,
        py::arg("cutoff")=DEFAULT_CUTOFF,
        py::arg("maxIter")=DEFAULT_MAX_ITER,
        py::arg("rtol")=DEFAULT_RTOL,
        py::arg("atol")=DEFAULT_ATOL);

    // Initialize all submodules
    init_XA(m);
    init_YA(m);
    init_YB(m);
    init_XC(m);
    init_YC(m);
    init_Minfinity(m);  // Add Minfinity initialization
}
