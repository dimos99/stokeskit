#ifndef XA_H
#define XA_H

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
 * @brief Calculates the pure number P_{n,p,q} for the XA scalar function in a two-body resistance problem.
 * 
 * @details
 * This function computes the P_{n,p,q} coefficients used in the series expansion of the two-body problem solutions.
 * (Eq. 3.6, 3.9 in Jeffrey & Onishi, 1981)
 * 
 * @param n The first index parameter
 * @param p The second index parameter
 * @param q The third index parameter
 * 
 * @return The value of the P_{n,p,q} number
 * 
 * @note Uses a global cache to store previously calculated values
 */
double P_npq_XA(int n, int p, int q);

/**
 * @brief Calculates the pure number V_{n,p,q} for the XA scalar function in a two-body resistance problem.
 * 
 * @details
 * This function computes the V_{n,p,q} coefficients used in the series expansion of the two-body problem solutions.
 * (Eq. 3.7, 3.8 in Jeffrey & Onishi, 1981)
 * 
 * @param n The first index parameter
 * @param p The second index parameter
 * @param q The third index parameter
 * 
 * @return The value of the V_{n,p,q} number
 * 
 * @note Uses a global cache to store previously calculated values
 */
double V_npq_XA(int n, int p, int q);


/**
 * @brief Calculates the function f_{k}(\lambda) for the XA scalar function in a two-body resistance problem.
 * 
 * @details
 * This function computes the f_{k}(\lambda) coefficients used in the series expansion of the two-body problem solutions.
 * (Eq. 3.15 in Jeffrey & Onishi, 1981)
 * 
 * @param k The order of the function
 * @param l The ratio of sphere radii
 * 
 * @return The value of the f_{k}(\lambda) function
 * 
 * @note Uses a global cache to store previously calculated values
 * @see P_npq_XA()
 */
double fk_XA(int k, double l);

/**
 * @brief Calculates the XA11 lubrication function for two spheres.
 * 
 * @details
 * This function computes the XA11 lubrication function which describes the hydrodynamic
 * interaction between two spheres. See Eq. 3.17 of Jeffrey and Onishi (1984) or Eq. 16 of Townsend (2023).
 * The normalization factor follows the one in Townsend (2023), or in Kim and Karrila (2005).
 * This is expected to be accurate for small separations only.
 * 
 * @param s Dimensionless surface-to-surface separation distance between spheres (s = 2*r/(a1+a2))
 * @param l Size ratio between spheres (l = a2/a1)
 * @param maxIter Maximum number of iterations for convergence in related calculations
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * 
 * @return The value of the XA11 lubrication function
 * 
 * @note The function uses cached g-coefficients for efficiency and combines them with
 *       logarithmic terms to compute the final result.
 */
double XA11_lubrication(double s, double l, int maxIter, double rtol, double atol);


/**
 * @brief Calculates the XA11 near-field function for two spheres.
 * 
 * @details
 * This function computes the XA11 lubrication function which describes the hydrodynamic
 * interaction between two spheres. See Eq. 3.20 of Jeffrey and Onishi (1984).
 * The normalization factor follows the one in Townsend (2023), or in Kim and Karrila (2005).
 * This is expected to be accurate for all separations, but is most efficient for small separations.
 * 
 * @param s Dimensionless surface-to-surface separation distance between spheres (s = 2*r/(a1+a2))
 * @param l Size ratio between spheres (l = a2/a1)
 * @param maxIter Maximum number of iterations for convergence in related calculations
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * 
 * @return The value of the XA11 lubrication function
 * 
 * @note The function uses cached g-coefficients for efficiency and combines them with
 *       logarithmic terms to compute the final result.
 */
double XA11_nf(double s, double l, int maxIter, double rtol, double atol);

/**
 * @brief Calculates the XA11 far-field function for two spheres.
 * 
 * @details
 * This function computes the XA11 lubrication function which describes the hydrodynamic
 * interaction between two spheres. See Eq. 3.13 of Jeffrey and Onishi (1984).
 * The normalization factor follows the one in Townsend (2023), or in Kim and Karrila (2005).
 * This is expected to be accurate for all separations, but is most efficient for small separations.
 * 
 * @param s Dimensionless surface-to-surface separation distance between spheres (s = 2*r/(a1+a2))
 * @param l Size ratio between spheres (l = a2/a1)
 * @param maxIter Maximum number of iterations for convergence in related calculations
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * 
 * @return The value of the XA11 lubrication function
 * 
 * @note The function uses cached g-coefficients for efficiency and combines them with
 *       logarithmic terms to compute the final result.
 */
double XA11_ff(double s, double l, int maxIter, double rtol, double atol);

/**
 * @brief Calculates the XA12 lubrication function for two spheres.
 * 
 * @details
 * This function computes the XA12 lubrication function which describes the hydrodynamic
 * interaction between two spheres. See Eq. 3.18 of Jeffrey and Onishi (1984) or Eq. 17 of Townsend (2023).
 * The normalization factor follows the one in Townsend (2023), or in Kim and Karrila (2005).
 * This is expected to be accurate for small separations only.
 * 
 * @param s Dimensionless surface-to-surface separation distance between spheres (s = 2*r/(a1+a2))
 * @param l Size ratio between spheres (l = a2/a1)
 * @param maxIter Maximum number of iterations for convergence in related calculations
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * 
 * @return The value of the XA12 lubrication function
 * 
 * @note The function uses cached g-coefficients for efficiency and combines them with
 *       logarithmic terms to compute the final result.
 */
double XA12_lubrication(double s, double l, int maxIter, double rtol, double atol);

/**
 * @brief Calculates the XA12 near-field function for two spheres.
 * 
 * @details
 * This function computes the XA12 lubrication function which describes the hydrodynamic
 * interaction between two spheres. See Eq. 3.21 of Jeffrey and Onishi (1984).
 * The normalization factor follows the one in Townsend (2023), or in Kim and Karrila (2005).
 * This is expected to be accurate for all separations, but is most efficient for small separations.
 * 
 * @param s Dimensionless surface-to-surface separation distance between spheres (s = 2*r/(a1+a2))
 * @param l Size ratio between spheres (l = a2/a1)
 * @param maxIter Maximum number of iterations for convergence in related calculations
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * 
 * @return The value of the XA12 lubrication function
 * 
 * @note The function uses cached g-coefficients for efficiency and combines them with
 *       logarithmic terms to compute the final result.
 */
double XA12_nf(double s, double l, int maxIter, double rtol, double atol);

/**
 * @brief Calculates the XA12 far-field function for two spheres.
 * 
 * @details
 * This function computes the XA12 lubrication function which describes the hydrodynamic
 * interaction between two spheres. See Eq. 3.14 of Jeffrey and Onishi (1984).
 * The normalization factor follows the one in Townsend (2023), or in Kim and Karrila (2005).
 * This is expected to be accurate for all separations, but is most efficient for small separations.
 * 
 * @param s Dimensionless surface-to-surface separation distance between spheres (s = 2*r/(a1+a2))
 * @param l Size ratio between spheres (l = a2/a1)
 * @param maxIter Maximum number of iterations for convergence in related calculations
 * @param rtol Relative tolerance for convergence
 * @param atol Absolute tolerance for convergence
 * 
 * @return The value of the XA12 lubrication function
 * 
 * @note The function uses cached g-coefficients for efficiency and combines them with
 *       logarithmic terms to compute the final result.
 */
double XA12_ff(double s, double l, int maxIter, double rtol, double atol);

/**
 * @brief Calculates the XA11 scalar resistance function
 * 
 * @details
 * This function computes the XA11 scalar resistance function which describes the hydrodynamic
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
 * @return double The XA11 function value
 */
double XA11(double s, double l, double lubr_cutoff, double cutoff, int maxIter, double rtol, double atol);

/**
 * @brief Calculates the XA12 scalar resistance function
 * 
 * @details
 * This function computes the XA12 scalar resistance function which describes the hydrodynamic
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
 * @return double The XA12 function value
 */
double XA12(double s, double l, double lubr_cutoff, double cutoff, int maxIter, double rtol, double atol);

// Add these declarations:
namespace xa_utils {
    void clear_cache();
}

#endif // XA_H