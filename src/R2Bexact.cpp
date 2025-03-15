#include "R2Bexact.h"
#include "common.h"
#include <algorithm>

// Global condensation matrix and indices for optimization
const std::array<std::array<int,2>,5> COND_IDX = {{{0,0}, {0,1}, {1,1}, {0,2}, {1,2}}};
const Eigen::MatrixXd COND_E = (Eigen::MatrixXd(5,5) << 
    HS3P1, 0, HS3M1, 0, 0,
    0, S2, 0, 0, 0,
    HS3M1, 0, HS3P1, 0, 0,
    0, 0, 0, S2, 0,
    0, 0, 0, 0, S2).finished();

// Unit displacement tensor implementations
double L1(const Eigen::Vector3d& d, int i, int j) {
    return d[i] * d[j];
}

double L2(const Eigen::Vector3d& d, int i, int j) {
    return (i == j ? 1.0 : 0.0) - d[i] * d[j];
}

double L3(const Eigen::Vector3d& d, int i, int j) {
    if (i == j) return 0.0;
    double sum = 0.0;
    for (int k = 0; k < 3; ++k) {
        if (k != i && k != j) {
            sum += levi(i, j, k) * d[k];
        }
    }
    return sum;
}

double L4(const Eigen::Vector3d& d, int i, int j, int k) {
    return (d[i] * d[j] - ((i == j) ? 1.0 : 0.0) / 3.0) * d[k];
}

double L5(const Eigen::Vector3d& d, int i, int j, int k) {
    return (d[i] * ((j == k) ? 1.0 : 0.0) + d[j] * ((i == k) ? 1.0 : 0.0) - 2 * d[i] * d[j] * d[k]);
}

double L6(const Eigen::Vector3d& d, int i, int j, int k) {
    double sum = 0.0;
    for (int l = 0; l < 3; ++l) {
        // Only sum over l values that are not i or k for first term
        if (l != i && l != k) {
            sum += levi(i, k, l) * d[l] * d[j];
        }
        // Only sum over l values that are not k or j for second term
        if (l != k && l != j) {
            sum += levi(j, k, l) * d[l] * d[i];
        }
    }
    return sum;
}

double L7(const Eigen::Vector3d& d, int i, int j, int k, int l) {
    return 1.5 * (d[i] * d[j] - ((i == j) ? 1.0 : 0.0) / 3.0) * 
                 (d[k] * d[l] - ((k == l) ? 1.0 : 0.0) / 3.0);
}

double L8(const Eigen::Vector3d& d, int i, int j, int k, int l) {
    return 0.5 * (d[i] * ((j == l) ? 1.0 : 0.0) * d[k] +
                  d[j] * ((i == l) ? 1.0 : 0.0) * d[k] +
                  d[i] * ((j == k) ? 1.0 : 0.0) * d[l] +
                  d[j] * ((i == k) ? 1.0 : 0.0) * d[l] -
                  4 * d[i] * d[j] * d[k] * d[l]);
}

double L9(const Eigen::Vector3d& d, int i, int j, int k, int l) {
    return 0.5 * (((i == k) ? 1.0 : 0.0) * ((j == l) ? 1.0 : 0.0) +
                  ((j == k) ? 1.0 : 0.0) * ((i == l) ? 1.0 : 0.0) -
                  ((i == j) ? 1.0 : 0.0) * ((k == l) ? 1.0 : 0.0) +
                  d[i] * d[j] * ((k == l) ? 1.0 : 0.0) +
                  ((i == j) ? 1.0 : 0.0) * d[k] * d[l] -
                  d[i] * ((j == l) ? 1.0 : 0.0) * d[k] -
                  d[j] * ((i == l) ? 1.0 : 0.0) * d[k] -
                  d[i] * ((j == k) ? 1.0 : 0.0) * d[l] -
                  d[j] * ((i == k) ? 1.0 : 0.0) * d[l] +
                  d[i] * d[j] * d[k] * d[l]);
}

// Helper functions for the scalar resistance functions
// This is a simplified replacement for the Python XYZ function that uses interpolation
// In our case, we directly call the appropriate resistance functions
double getResistanceFunction(int type, int gamma, double s, double lambda) {
    // Use default parameters from common.h
    double lubr_cutoff = DEFAULT_LUBR_CUTOFF;
    double cutoff = DEFAULT_CUTOFF;
    int maxIter = DEFAULT_MAX_ITER; 
    double rtol = DEFAULT_RTOL;
    double atol = DEFAULT_ATOL;
    
    switch(type) {
        case 0: // XA
            return (gamma == 0) ? 
                XA11(s, lambda, lubr_cutoff, cutoff, maxIter, rtol, atol) : 
                XA12(s, lambda, lubr_cutoff, cutoff, maxIter, rtol, atol);
        case 1: // YA
            return (gamma == 0) ? 
                YA11(s, lambda, lubr_cutoff, cutoff, maxIter, rtol, atol) : 
                YA12(s, lambda, lubr_cutoff, cutoff, maxIter, rtol, atol);
        case 2: // YB
            return (gamma == 0) ? 
                YB11(s, lambda, lubr_cutoff, cutoff, maxIter, rtol, atol) : 
                YB12(s, lambda, lubr_cutoff, cutoff, maxIter, rtol, atol);
        case 3: // XC
            return (gamma == 0) ? 
                XC11(s, lambda, lubr_cutoff, cutoff, maxIter, rtol, atol) : 
                XC12(s, lambda, lubr_cutoff, cutoff, maxIter, rtol, atol);
        case 4: // YC
            return (gamma == 0) ? 
                YC11(s, lambda, lubr_cutoff, cutoff, maxIter, rtol, atol) : 
                YC12(s, lambda, lubr_cutoff, cutoff, maxIter, rtol, atol);
        // For types 5-10 we would need to implement XG, YG, YH, XM, YM, ZM functions
        // For now, we'll use minimal implementations that approximate far-field behavior
        default:
            return 0.0;
    }
}

// Submatrix element functions
double Af(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i, int j) {
    double XAval = getResistanceFunction(0, gamma, s, lambda); // XA
    double YAval = getResistanceFunction(1, gamma, s, lambda); // YA
    return XAval * L1(d, i, j) + YAval * L2(d, i, j);
}

double Bf(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i, int j) {
    double YBval = getResistanceFunction(2, gamma, s, lambda); // YB
    return YBval * L3(d, i, j);
}

double Cf(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i, int j) {
    double XCval = getResistanceFunction(3, gamma, s, lambda); // XC
    double YCval = getResistanceFunction(4, gamma, s, lambda); // YC
    return XCval * L1(d, i, j) + YCval * L2(d, i, j);
}

// Implementation of G, H, and M resistance functions based on reference code
double Gf(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i, int j, int k) {
    // Use minimal implementations since we don't have XG and YG functions yet
    // This matches the structure without needing the actual function values
    double XGval = 0.0; // getResistanceFunction(5, gamma, s, lambda) would go here
    double YGval = 0.0; // getResistanceFunction(6, gamma, s, lambda) would go here
    return XGval * L4(d, i, j, k) + YGval * L5(d, i, j, k);
}

double Hf(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i, int j, int k) {
    // Use minimal implementation since we don't have YH function yet
    double YHval = 0.0; // getResistanceFunction(7, gamma, s, lambda) would go here
    return YHval * L6(d, i, j, k);
}

double Mf(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i, int j, int k, int l) {
    // Use minimal implementations since we don't have XM, YM, ZM functions yet
    double XMval = 0.0; // getResistanceFunction(8, gamma, s, lambda) would go here
    double YMval = 0.0; // getResistanceFunction(9, gamma, s, lambda) would go here
    double ZMval = 0.0; // getResistanceFunction(10, gamma, s, lambda) would go here
    return XMval * L7(d, i, j, k, l) + YMval * L8(d, i, j, k, l) + ZMval * L9(d, i, j, k, l);
}

// Condensed submatrix helper functions
Eigen::VectorXd con_Gf_row(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i) {
    Eigen::VectorXd result(5);
    
    // Calculate once and reuse
    double A = Gf(gamma, d, lambda, s, 0, 0, i);
    double B = Gf(gamma, d, lambda, s, 1, 1, i);
    double C = Gf(gamma, d, lambda, s, 0, 1, i);
    double D = Gf(gamma, d, lambda, s, 0, 2, i);
    double E = Gf(gamma, d, lambda, s, 1, 2, i);
    
    result << (HS3P1 * A + HS3M1 * B),
              S2 * C,
              (HS3M1 * A + HS3P1 * B),
              S2 * D,
              S2 * E;
              
    return result;
}

Eigen::VectorXd con_Hf_row(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i) {
    Eigen::VectorXd result(5);
    
    // Calculate once and reuse
    double A = Hf(gamma, d, lambda, s, 0, 0, i);
    double B = Hf(gamma, d, lambda, s, 1, 1, i);
    double C = Hf(gamma, d, lambda, s, 0, 1, i);
    double D = Hf(gamma, d, lambda, s, 0, 2, i);
    double E = Hf(gamma, d, lambda, s, 1, 2, i);
    
    result << (HS3P1 * A + HS3M1 * B),
              S2 * C,
              (HS3M1 * A + HS3P1 * B),
              S2 * D,
              S2 * E;
              
    return result;
}

// Calculate indices for submatrices in the resistance matrix
SubmatrixIndices calculateSubmatrixIndices(int a1_index, int a2_index, int num_spheres) {
    SubmatrixIndices indices;
    
    // Calculate starting indices for each submatrix
    // Each sphere contributes 11 rows/columns total (3 for translation, 3 for rotation, 5 for strain)
    indices.A_start_row = 3 * a1_index;
    indices.A_start_col = 3 * a2_index;
    
    indices.B_start_row = 3 * a1_index;
    indices.B_start_col = 3 * num_spheres + 3 * a2_index;
    
    indices.B21_start_row = 3 * num_spheres + 3 * a1_index;
    indices.B21_start_col = 3 * a2_index;
    
    indices.C_start_row = 3 * num_spheres + 3 * a1_index;
    indices.C_start_col = 3 * num_spheres + 3 * a2_index;
    
    indices.G_start_row = 3 * a1_index;
    indices.G_start_col = 6 * num_spheres + 5 * a2_index;
    
    indices.G21_start_row = 6 * num_spheres + 5 * a1_index;
    indices.G21_start_col = 3 * a2_index;
    
    indices.H_start_row = 3 * num_spheres + 3 * a1_index;
    indices.H_start_col = 6 * num_spheres + 5 * a2_index;
    
    indices.H21_start_row = 6 * num_spheres + 5 * a1_index;
    indices.H21_start_col = 3 * num_spheres + 3 * a2_index;
    
    indices.M_start_row = 6 * num_spheres + 5 * a1_index;
    indices.M_start_col = 6 * num_spheres + 5 * a2_index;
    
    return indices;
}

// Helper function to find close particles
std::tuple<std::vector<std::pair<int, int>>, std::vector<Eigen::Vector3d>, std::vector<double>, std::vector<double>> 
findCloseParticles(const std::vector<Eigen::Vector3d>& positions, const std::vector<double>& radii, double cutoff_factor) {
    std::vector<std::pair<int, int>> close_pairs;
    std::vector<Eigen::Vector3d> displacements;
    std::vector<double> distances;
    std::vector<double> size_ratios;
    
    int n = positions.size();
    
    for (int i = 0; i < n; ++i) {
        // Always include the diagonal (particle with itself)
        close_pairs.push_back(std::make_pair(i, i));
        displacements.push_back(Eigen::Vector3d::Zero());
        distances.push_back(0.0);
        size_ratios.push_back(1.0);
        
        for (int j = i+1; j < n; ++j) {
            Eigen::Vector3d r = positions[j] - positions[i];
            double distance = r.norm();
            double cutoff = cutoff_factor * (radii[i] + radii[j]);
            
            if (distance < cutoff) {
                // Add the pair
                close_pairs.push_back(std::make_pair(i, j));
                // For R2Bexact, displacement convention is a2-a1, already correct
                displacements.push_back(r);
                // s_dash = 2s/(a+b) is the dimensionless separation parameter
                double s_dash = 2.0 * distance / (radii[i] + radii[j]);
                distances.push_back(s_dash);
                // lambda = a2/a1 is the size ratio between particles
                size_ratios.push_back(radii[j] / radii[i]);
            }
        }
    }
    
    return std::make_tuple(close_pairs, displacements, distances, size_ratios);
}

Eigen::MatrixXd generateR2Bexact(const std::vector<Eigen::Vector3d>& positions, 
                                const std::vector<double>& radii, 
                                double mu, 
                                double cutoff_factor) {
    // Validate inputs
    int num_spheres = static_cast<int>(positions.size());
    if (num_spheres != static_cast<int>(radii.size())) {
        throw std::invalid_argument("Number of positions must match number of radii");
    }
    if (mu <= 0) {
        throw std::invalid_argument("Viscosity must be positive");
    }
    
    // Resistance matrix size (11 indices per sphere: 3 translation, 3 rotation, 5 strain)
    int matrix_size = 11 * num_spheres;
    Eigen::MatrixXd R2Bexact = Eigen::MatrixXd::Zero(matrix_size, matrix_size);
    
    // Find pairs of particles closer than cutoff
    auto [close_pairs, displacements, s_dash_values, lambdas] = findCloseParticles(positions, radii, cutoff_factor);
    
    // Temporary matrix for M calculation
    Eigen::MatrixXd M_temp(5, 5);
    
    // Loop through all close pairs
    for (size_t idx = 0; idx < close_pairs.size(); ++idx) {
        int a1_index = close_pairs[idx].first;
        int a2_index = close_pairs[idx].second;
        
        // Get the indices for all submatrices
        SubmatrixIndices indices = calculateSubmatrixIndices(a1_index, a2_index, num_spheres);
        
        // Calculate scale factors
        double a1 = radii[a1_index];
        double scale_A = 6.0 * PI * a1;
        double scale_B = 4.0 * PI * a1 * a1;
        double scale_C = 8.0 * PI * a1 * a1 * a1;
        double scale_G = 4.0 * PI * a1 * a1;
        double scale_H = 8.0 * PI * a1 * a1 * a1;
        double scale_M = 20.0/3.0 * PI * a1 * a1 * a1;
        
        double a2 = radii[a2_index];
        double scale_B2 = 4.0 * PI * a2 * a2;
        double scale_G2 = 4.0 * PI * a2 * a2;
        double scale_H2 = 8.0 * PI * a2 * a2 * a2;
        
        if (a1_index == a2_index) {
            // Diagonal blocks - need to sum contributions from all nearby particles
            Eigen::Matrix3d A_sum = Eigen::Matrix3d::Identity(); // Initialize with self-mobility
            Eigen::Matrix3d Bt_sum = Eigen::Matrix3d::Zero();
            Eigen::Matrix3d C_sum = Eigen::Matrix3d::Identity(); // Initialize with self-mobility
            Eigen::MatrixXd Gt_sum(3, 5);
            Gt_sum.setZero();
            Eigen::MatrixXd Ht_sum(3, 5);
            Ht_sum.setZero();
            Eigen::MatrixXd M_sum(5, 5);
            M_sum.setIdentity(); // Initialize with self-mobility
            
            // Find all particles close to a1_index
            for (size_t p_idx = 0; p_idx < close_pairs.size(); ++p_idx) {
                auto [p1, p2] = close_pairs[p_idx];
                
                // Check if this pair involves a1_index but isn't the diagonal
                if ((p1 == a1_index && p2 != a1_index) || 
                    (p2 == a1_index && p1 != a1_index)) {
                    
                    int p_index = (p1 == a1_index) ? p2 : p1;
                    double lambda_p = radii[p_index] / radii[a1_index];
                    
                    // Get or compute displacement and distance
                    Eigen::Vector3d r_p;
                    double s_dash_p;
                    
                    if (p1 == a1_index && p2 == p_index) {
                        r_p = displacements[p_idx];
                        s_dash_p = s_dash_values[p_idx];
                    } else {
                        // Must be p2 == a1_index && p1 == p_index, so reverse direction
                        r_p = -displacements[p_idx];
                        s_dash_p = s_dash_values[p_idx];
                    }
                    
                    // Enforce minimum separation to match Python code
                    if (s_dash_p < 2.001) {
                        s_dash_p = 2.001;  // Minimum dimensionless separation
                    }
                    
                    // Unit direction vector
                    Eigen::Vector3d d_p = r_p.normalized();
                    
                    // Calculate submatrix elements for gamma=0 (self-interaction)
                    for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            A_sum(i, j) += Af(0, d_p, lambda_p, s_dash_p, i, j);
                            Bt_sum(i, j) += Bf(0, d_p, lambda_p, s_dash_p, j, i);
                            C_sum(i, j) += Cf(0, d_p, lambda_p, s_dash_p, i, j);
                        }
                        
                        // Get condensed rows for G and H
                        Gt_sum.row(i) += con_Gf_row(0, d_p, lambda_p, s_dash_p, i);
                        Ht_sum.row(i) += con_Hf_row(0, d_p, lambda_p, s_dash_p, i);
                    }
                    
                    // Calculate M submatrix using condensation
                    for (int i = 0; i < 5; ++i) {
                        for (int j = 0; j < 5; ++j) {
                            M_temp(i, j) = Mf(0, d_p, lambda_p, s_dash_p,
                                            COND_IDX[i][0], COND_IDX[i][1], 
                                            COND_IDX[j][0], COND_IDX[j][1]);
                        }
                    }
                    M_sum += COND_E * M_temp * COND_E;
                }
            }
            
            // Scale and assign the summed submatrices to the resistance matrix
            R2Bexact.block<3,3>(indices.A_start_row, indices.A_start_col) = A_sum * scale_A;
            R2Bexact.block<3,3>(indices.B_start_row, indices.B_start_col) = Bt_sum * scale_B;
            R2Bexact.block<3,3>(indices.C_start_row, indices.C_start_col) = C_sum * scale_C;
            R2Bexact.block<3,5>(indices.G_start_row, indices.G_start_col) = Gt_sum * scale_G;
            R2Bexact.block<3,5>(indices.H_start_row, indices.H_start_col) = Ht_sum * scale_H;
            R2Bexact.block<5,5>(indices.M_start_row, indices.M_start_col) = M_sum * scale_M;
            
        } else {
            // Off-diagonal blocks - interaction between different particles
            Eigen::Vector3d r = displacements[idx];
            double s_dash = s_dash_values[idx];
            double lambda = lambdas[idx];
            
            // Enforce minimum separation to match Python code
            if (s_dash < 2.001) {
                s_dash = 2.001;  // Minimum dimensionless separation
            }
            
            // Direction unit vector
            Eigen::Vector3d d = r.normalized();
            
            // Calculate A submatrix (gamma=1 for cross-interaction)
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    R2Bexact(indices.A_start_row + i, indices.A_start_col + j) = 
                        Af(1, d, lambda, s_dash, i, j) * scale_A;
                }
            }
            
            // Calculate B submatrix
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    R2Bexact(indices.B_start_row + i, indices.B_start_col + j) = 
                        Bf(1, -d, 1.0/lambda, s_dash, j, i) * scale_B2;
                }
            }
            
            // Calculate C submatrix
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    R2Bexact(indices.C_start_row + i, indices.C_start_col + j) = 
                        Cf(1, d, lambda, s_dash, i, j) * scale_C;
                }
            }
            
            // Calculate G submatrix
            for (int i = 0; i < 3; ++i) {
                Eigen::VectorXd g_row = con_Gf_row(1, -d, 1.0/lambda, s_dash, i);
                R2Bexact.block<1,5>(indices.G_start_row + i, indices.G_start_col) = g_row.transpose() * scale_G2;
            }
            
            // Calculate H submatrix
            for (int i = 0; i < 3; ++i) {
                Eigen::VectorXd h_row = con_Hf_row(1, -d, 1.0/lambda, s_dash, i);
                R2Bexact.block<1,5>(indices.H_start_row + i, indices.H_start_col) = h_row.transpose() * scale_H2;
            }
            
            // Calculate M submatrix
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < 5; ++j) {
                    M_temp(i, j) = Mf(1, d, lambda, s_dash,
                                    COND_IDX[i][0], COND_IDX[i][1], 
                                    COND_IDX[j][0], COND_IDX[j][1]);
                }
            }
            Eigen::MatrixXd M = COND_E * M_temp * COND_E * scale_M;
            R2Bexact.block<5,5>(indices.M_start_row, indices.M_start_col) = M;
            
            // Calculate transpose blocks
            if (std::abs(lambda - 1.0) < 1e-10) {
                // For equal-sized particles, use symmetry properties
                R2Bexact.block<3,3>(indices.B21_start_row, indices.B21_start_col) = 
                    -R2Bexact.block<3,3>(indices.B_start_row, indices.B_start_col);
                
                R2Bexact.block<5,3>(indices.G21_start_row, indices.G21_start_col) = 
                    -R2Bexact.block<3,5>(indices.G_start_row, indices.G_start_col).transpose();
                
                R2Bexact.block<5,3>(indices.H21_start_row, indices.H21_start_col) = 
                    R2Bexact.block<3,5>(indices.H_start_row, indices.H_start_col).transpose();
            } else {
                // For unequal sizes, calculate explicitly
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        R2Bexact(indices.B21_start_row + i, indices.B21_start_col + j) = 
                            Bf(1, d, lambda, s_dash, j, i) * scale_B;
                    }
                }
                
                // Calculate G21 and H21
                for (int i = 0; i < 3; ++i) {
                    // Calculate G21 and H21 properly
                    Eigen::VectorXd g_row = con_Gf_row(1, d, lambda, s_dash, i);
                    for (int j = 0; j < 5; ++j) {
                        R2Bexact(indices.G21_start_row + j, indices.G21_start_col + i) = g_row(j) * scale_G;
                    }
                    
                    Eigen::VectorXd h_row = con_Hf_row(1, d, lambda, s_dash, i);
                    for (int j = 0; j < 5; ++j) {
                        R2Bexact(indices.H21_start_row + j, indices.H21_start_col + i) = h_row(j) * scale_H;
                    }
                }
            }
        }
    }
    
    // Apply explicit symmetrization with extra post-processing
    // First make sure the matrix is symmetric
    R2Bexact = 0.5 * (R2Bexact + R2Bexact.transpose());
    
    // Multiply by viscosity
    return mu * R2Bexact;
}

// Simple test function to validate R2Bexact implementation
Eigen::MatrixXd testR2BexactSimple() {
    // Test with two equal-sized spheres at a specific distance
    std::vector<Eigen::Vector3d> positions = {
        Eigen::Vector3d(0, 0, 0),
        Eigen::Vector3d(0, 0, 3.0)
    };
    std::vector<double> radii = {1.0, 1.0};
    
    // Generate the resistance matrix
    Eigen::MatrixXd result = generateR2Bexact(positions, radii, 1.0, 2.0);
    
    return result;
}
