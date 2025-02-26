#ifndef YC_H
#define YC_H

#pragma once

#include <cmath>       // For std::pow, std::log, std::abs
#include <iostream>    // For std::cout
#include <iomanip>     // For std::setw
#include <map>         // For caching
#include <chrono>
#include <unordered_map>
#include <string>
#include <mutex>
#include "common.h"

// Extra includes
#include <pybind11/pybind11.h>  // Main pybind11 header
#include <pybind11/stl.h>       // For std::vector bindings
#include <pybind11/numpy.h>     // For numpy bindings


namespace py = pybind11;



// Declarations for functions

/**
 * @brief Calculates the pure number P_{n,p,q} for the YC scalar function in a two-body resistance problem.
 * 
 * @details
 * This function computes the P_{n,p,q} coefficients used in the series expansion of the two-body problem solutions.
 * (See Jeffrey & Onishi, 1981)
 * 
 * @param n The first index parameter
 * @param p The second index parameter
 * @param q The third index parameter
 * 
 * @return The value of the P_{n,p,q} number
 * 
 * @note Uses a global cache to store previously calculated values
 */
double P_npq_YC(int n, int p, int q);

/**
 * @brief Calculates the pure number V_{n,p,q} for the YC scalar function in a two-body resistance problem.
 * 
 * @details
 * This function computes the V_{n,p,q} coefficients used in the series expansion of the two-body problem solutions.
 * (See Jeffrey & Onishi, 1981)
 * 
 * @param n The first index parameter
 * @param p The second index parameter
 * @param q The third index parameter
 * 
 * @return The value of the V_{n,p,q} number
 * 
 * @note Uses a global cache to store previously calculated values
 */
double V_npq_YC(int n, int p, int q);

/**
 * @brief Calculates the pure number Q_{n,p,q} for the YC scalar function in a two-body resistance problem.
 * 
 * @details
 * This function computes the Q_{n,p,q} coefficients used in the series expansion of the two-body problem solutions.
 * (See Jeffrey & Onishi, 1981)
 * 
 * @param n The first index parameter
 * @param p The second index parameter
 * @param q The third index parameter
 * 
 * @return The value of the Q_{n,p,q} number
 * 
 * @note Uses a global cache to store previously calculated values
 */
double Q_npq_YC(int n, int p, int q);


/**
 * @brief Calculates the function f_{k}(\lambda) for the YC scalar function in a two-body resistance problem.
 * 
 * @details
 * This function computes the f_{k}(\lambda) coefficients used in the series expansion of the two-body problem solutions.
 * (See Jeffrey & Onishi, 1981)
 * 
 * @param k The order of the function
 * @param l The ratio of sphere radii
 * 
 * @return The value of the f_{k}(\lambda) function
 * 
 * @note Uses a global cache to store previously calculated values
 * @see P_npq_YC()
 */
double fk_YC(int k, double l);

/**
 * @brief Calculates the YC11 lubrication function for two spheres.
 * 
 * @details
 * This function computes the YC11 lubrication function which describes the hydrodynamic
 * interaction between two spheres. See Jeffrey and Onishi (1984) or Townsend (2023).
 * The normalization factor follows the one in Townsend (2023), or in Kim and Karrila (2005).
 * This is expected to be accurate for small separations only.
 * 
 * @param s Dimensionless surface-to-surface separation distance between spheres (s = 2*r/(a1+a2))
 * @param l Size ratio between spheres (l = a2/a1)
 * @param maxIter Maximum number of iterations for convergence in related calculations
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * 
 * @return The value of the YC11 lubrication function
 * 
 * @note The function uses cached g-coefficients for efficiency and combines them with
 *       logarithmic terms to compute the final result.
 */
double YC11_lubrication(double s, double l, int maxIter, double rtol, double atol);


/**
 * @brief Calculates the YC11 near-field function for two spheres.
 * 
 * @details
 * This function computes the YC11 lubrication function which describes the hydrodynamic
 * interaction between two spheres. See Jeffrey and Onishi (1984).
 * The normalization factor follows the one in Townsend (2023), or in Kim and Karrila (2005).
 * This is expected to be accurate for all separations, but is most efficient for small separations.
 * 
 * @param s Dimensionless surface-to-surface separation distance between spheres (s = 2*r/(a1+a2))
 * @param l Size ratio between spheres (l = a2/a1)
 * @param maxIter Maximum number of iterations for convergence in related calculations
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * 
 * @return The value of the YC11 lubrication function
 * 
 * @note The function uses cached g-coefficients for efficiency and combines them with
 *       logarithmic terms to compute the final result.
 */
double YC11_nf(double s, double l, int maxIter, double rtol, double atol);

/**
 * @brief Calculates the YC11 far-field function for two spheres.
 * 
 * @details
 * This function computes the YC11 lubrication function which describes the hydrodynamic
 * interaction between two spheres. See Jeffrey and Onishi (1984).
 * The normalization factor follows the one in Townsend (2023), or in Kim and Karrila (2005).
 * This is expected to be accurate for all separations, but is most efficient for small separations.
 * 
 * @param s Dimensionless surface-to-surface separation distance between spheres (s = 2*r/(a1+a2))
 * @param l Size ratio between spheres (l = a2/a1)
 * @param maxIter Maximum number of iterations for convergence in related calculations
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * 
 * @return The value of the YC11 lubrication function
 * 
 * @note The function uses cached g-coefficients for efficiency and combines them with
 *       logarithmic terms to compute the final result.
 */
double YC11_ff(double s, double l, int maxIter, double rtol, double atol);

/**
 * @brief Calculates the YC12 lubrication function for two spheres.
 * 
 * @details
 * This function computes the YC12 lubrication function which describes the hydrodynamic
 * interaction between two spheres. See Jeffrey and Onishi (1984) or Townsend (2023).
 * The normalization factor follows the one in Townsend (2023), or in Kim and Karrila (2005).
 * This is expected to be accurate for small separations only.
 * 
 * @param s Dimensionless surface-to-surface separation distance between spheres (s = 2*r/(a1+a2))
 * @param l Size ratio between spheres (l = a2/a1)
 * @param maxIter Maximum number of iterations for convergence in related calculations
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * 
 * @return The value of the YC12 lubrication function
 * 
 * @note The function uses cached g-coefficients for efficiency and combines them with
 *       logarithmic terms to compute the final result.
 */
double YC12_lubrication(double s, double l, int maxIter, double rtol, double atol);

/**
 * @brief Calculates the YC12 near-field function for two spheres.
 * 
 * @details
 * This function computes the YC12 lubrication function which describes the hydrodynamic
 * interaction between two spheres. See Jeffrey and Onishi (1984).
 * The normalization factor follows the one in Townsend (2023), or in Kim and Karrila (2005).
 * This is expected to be accurate for all separations, but is most efficient for small separations.
 * 
 * @param s Dimensionless surface-to-surface separation distance between spheres (s = 2*r/(a1+a2))
 * @param l Size ratio between spheres (l = a2/a1)
 * @param maxIter Maximum number of iterations for convergence in related calculations
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * 
 * @return The value of the YC12 lubrication function
 * 
 * @note The function uses cached g-coefficients for efficiency and combines them with
 *       logarithmic terms to compute the final result.
 */
double YC12_nf(double s, double l, int maxIter, double rtol, double atol);

/**
 * @brief Calculates the YC12 far-field function for two spheres.
 * 
 * @details
 * This function computes the YC12 lubrication function which describes the hydrodynamic
 * interaction between two spheres. See Jeffrey and Onishi (1984).
 * The normalization factor follows the one in Townsend (2023), or in Kim and Karrila (2005).
 * This is expected to be accurate for all separations, but is most efficient for small separations.
 * 
 * @param s Dimensionless surface-to-surface separation distance between spheres (s = 2*r/(a1+a2))
 * @param l Size ratio between spheres (l = a2/a1)
 * @param maxIter Maximum number of iterations for convergence in related calculations
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * 
 * @return The value of the YC12 lubrication function
 * 
 * @note The function uses cached g-coefficients for efficiency and combines them with
 *       logarithmic terms to compute the final result.
 */
double YC12_ff(double s, double l, int maxIter, double rtol, double atol);

/**
 * @brief Calculates the YC11 scalar resistance function
 * 
 * @details
 * This function computes the YC11 scalar resistance function which describes the hydrodynamic
 * interaction between two spheres. The function uses a combination of lubrication, near-field,
 * and far-field approximations to compute the final result. The function automatically selects
 * the appropriate approximation based on the input parameters.
 *
 * @param[in] s Dimensionless surface-to-surface separation
 * @param[in] l Ratio of sphere radii (a2/a1)
 * @param[in] lubr_cutoff Lubrication cutoff distance
 * @param[in] cutoff Far-field cutoff distance
 * @param[in] maxIter Maximum iterations for convergence
 * @param[in] rtol Relative tolerance
 * @param[in] atol Absolute tolerance
 * @return double The YC11 function value
 */
double YC11(double s, double l, double lubr_cutoff, double cutoff, int maxIter, double rtol, double atol);

/**
 * @brief Calculates the YC12 scalar resistance function
 * 
 * @details
 * This function computes the YC12 scalar resistance function which describes the hydrodynamic
 * interaction between two spheres. The function uses a combination of lubrication, near-field,
 * and far-field approximations to compute the final result. The function automatically selects
 * the appropriate approximation based on the input parameters.
 *
 * @param[in] s Dimensionless surface-to-surface separation
 * @param[in] l Ratio of sphere radii (a2/a1)
 * @param[in] lubr_cutoff Lubrication cutoff distance
 * @param[in] cutoff Far-field cutoff distance
 * @param[in] maxIter Maximum iterations for convergence
 * @param[in] rtol Relative tolerance
 * @param[in] atol Absolute tolerance
 * @return double The YC12 function value
 */
double YC12(double s, double l, double lubr_cutoff, double cutoff, int maxIter, double rtol, double atol);

// Add these declarations:
namespace yc_utils {
    void clear_cache();
}

#endif // YC_H