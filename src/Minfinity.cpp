#include "Minfinity.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

const double PI = 3.141592653589793;
const double S2 = std::sqrt(2.0);
const double S3 = std::sqrt(3.0);
const double HS3P1 = 0.5*(S3+1);
const double HS3M1 = 0.5*(S3-1);

using namespace Eigen;

// Constants for the COND_IDX and COND_E matrices
const std::array<std::array<int,2>,5> COND_IDX = {{{0,0}, {0,1}, {1,1}, {0,2}, {1,2}}};
const MatrixXd COND_E = (MatrixXd(5,5) << 
    HS3P1, 0, HS3M1, 0, 0,
    0, S2, 0, 0, 0,
    HS3M1, 0, HS3P1, 0, 0,
    0, 0, 0, S2, 0,
    0, 0, 0, 0, S2).finished();

// Helper functions for tensor operations
int levi(int i, int j, int k) {
    if (i == j || j == k || k == i) return 0;
    if ((i == 0 && j == 1 && k == 2) || 
        (i == 1 && j == 2 && k == 0) || 
        (i == 2 && j == 0 && k == 1)) return 1;
    return -1;
}

// J tensor
Matrix3d J_tensor(const Vector3d& r, double s) {
    Matrix3d J = Matrix3d::Zero();
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            J(i,j) = ((i == j ? 1.0 : 0.0) / s + r[i]*r[j]/(s*s*s));
        }
    }
    return J;
}

// Laplacian of J
Matrix3d Lap_J(const Vector3d& r, double s) {
    Matrix3d result = Matrix3d::Zero();
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            result(i,j) = 2.0 * (i == j ? 1.0 : 0.0)/(s*s*s) - 6.0*r[i]*r[j]/(s*s*s*s*s);
        }
    }
    return result;
}

// Fix D_J to handle multiple indices
Matrix3d D_J(const Vector3d& r, double s, int l, int i, int j) {
    Matrix3d result = Matrix3d::Zero();
    for(int m = 0; m < 3; m++) {
        for(int n = 0; n < 3; n++) {
            result(m,n) = (-((m == n ? 1.0 : 0.0)*r[l] - (m == l ? 1.0 : 0.0)*r[n] 
                          - (n == l ? 1.0 : 0.0)*r[m])/s/s/s 
                          + 3*r[m]*r[n]*r[l]/pow(s,5));
        }
    }
    return result;
}

// R tensor
Matrix3d R_tensor(const Vector3d& r, double s) {
    Matrix3d R = Matrix3d::Zero();
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            double sum = 0;
            for(int k = 0; k < 3; k++) {
                for(int l = 0; l < 3; l++) {
                    if(k != l && j != k && j != l) {
                        Matrix3d DJ = D_J(r,s,k,l,i); 
                        sum += levi(j,k,l) * DJ(0,0);
                    }
                }
            }
            R(i,j) = -0.5 * sum;
        }
    }
    return R;
}

// Fix K_tensor to properly handle matrix operations
Matrix3d K_tensor(const Vector3d& r, double s, int i, int j, int k) {
    Matrix3d K = Matrix3d::Zero();
    for(int m = 0; m < 3; m++) {
        for(int n = 0; n < 3; n++) {
            double v1 = D_J(r,s,k,n,m)(0,0); // Extract scalar value
            double v2 = D_J(r,s,n,m,k)(0,0); // Extract scalar value
            K(m,n) = 0.5 * (v1 + v2);
        }
    }
    return K;
}

// Fix DD_J to compute elements directly
Matrix3d DD_J(const Vector3d& r, double s, int m, int l, int i, int j) {
    Matrix3d result = Matrix3d::Zero();
    double s3 = pow(s,3);
    double s5 = pow(s,5); 
    double s7 = pow(s,7);

    // Compute each element directly
    for(int p = 0; p < 3; p++) {
        for(int q = 0; q < 3; q++) {
            result(p,q) = (-(p == q ? 1.0 : 0.0)*(l == m ? 1.0 : 0.0) + 
                           (p == l ? 1.0 : 0.0)*(q == m ? 1.0 : 0.0) + 
                           (p == m ? 1.0 : 0.0)*(q == l ? 1.0 : 0.0))/s3
                         - 3.0*(-(p == q ? 1.0 : 0.0)*r[l]*r[m] + 
                               (p == l ? 1.0 : 0.0)*r[q]*r[m] +
                               (p == m ? 1.0 : 0.0)*r[q]*r[l])/s5
                         + 15.0*r[p]*r[q]*r[l]*r[m]/s7;
        }
    }
    return result;
}

// Fix DLap_J calculation
Matrix3d DLap_J(const Vector3d& r, double s, int k, int i, int j) {
    double s5 = pow(s,5);
    double s7 = pow(s,7);
    
    Matrix3d result = Matrix3d::Zero();
    
    // First term 
    result = -6.0/s5 * r[k] * Matrix3d::Identity(); 
    
    // Second term - construct proper outer products
    Vector3d ei = Vector3d::Unit(i);
    Vector3d ej = Vector3d::Unit(j);
    result += -6.0/s5 * (ei * r[j] * ej.transpose() + ej * r[i] * ei.transpose());
    
    // Third term
    result += 30.0/s7 * r[i]*r[j]*r[k] * Matrix3d::Identity();
    
    return result;
}

Matrix3d Lap_R(const Vector3d& r, double s) {
    Matrix3d result = Matrix3d::Zero();
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            for(int k = 0; k < 3; k++) {
                for(int l = 0; l < 3; l++) {
                    if(k != l && j != k && j != l) {
                        Matrix3d DL = DLap_J(r,s,k,i,l);
                        result(i,j) += levi(j,k,l) * DL.trace();
                    }
                }
            }
        }
    }
    return -0.5 * result;
}

Matrix3d Lap_K(const Vector3d& r, double s, int i, int j, int k) {
    return DLap_J(r,s,i,j,k);
}

Matrix3d DLap_K(const Vector3d& r, double s, int l, int i, int j, int k) {
    double s5 = pow(s,5);
    double s7 = pow(s,7);
    double s9 = pow(s,9);
    
    Matrix3d result = Matrix3d::Zero();
    Vector3d term1 = Vector3d::Zero();
    
    // First term
    result = -6.0/s5 * (DD_J(r,s,i,j,k,l) + DD_J(r,s,i,k,j,l));
    
    // Second term
    result = result - 210.0/s9 * r[i]*r[j]*r[k]*r[l] * Matrix3d::Identity();
    
    // Third term
    double term3 = 30.0/s7 * (
        (i == j ? 1.0 : 0.0)*r[k]*r[l] +
        (i == k ? 1.0 : 0.0)*r[j]*r[l] +
        (j == k ? 1.0 : 0.0)*r[i]*r[l] +
        (i == l ? 1.0 : 0.0)*r[j]*r[k] +
        (j == l ? 1.0 : 0.0)*r[i]*r[k] +
        (k == l ? 1.0 : 0.0)*r[i]*r[j]
    );
    result = result + term3 * Matrix3d::Identity();
    
    return result;
}

// Add the matrix contraction helpers
Vector5d con_M13_row(const Vector3d& r, double s, double a1, double a2, int i, double c, double mu) {
    Vector5d result;
    if(s > 1e-10) {
        Matrix3d K00 = K_tensor(r,s,i,0,0);
        Matrix3d K11 = K_tensor(r,s,i,1,1);
        Matrix3d LK00 = Lap_K(r,s,i,0,0);
        Matrix3d LK11 = Lap_K(r,s,i,1,1);
        
        double A = -c*(K00(0,0) + (a1*a1/6.0 + a2*a2/10.0)*LK00(0,0));
        double B = -c*(K11(0,0) + (a1*a1/6.0 + a2*a2/10.0)*LK11(0,0));
        
        Matrix3d K01 = K_tensor(r,s,i,0,1);
        Matrix3d K02 = K_tensor(r,s,i,0,2);
        Matrix3d K12 = K_tensor(r,s,i,1,2);
        
        result << (HS3P1*A + HS3M1*B),
                  S2*(-c*(K01(0,0) + (a1*a1/6.0 + a2*a2/10.0)*Lap_K(r,s,i,0,1)(0,0))),
                  (HS3M1*A + HS3P1*B),
                  S2*(-c*(K02(0,0) + (a1*a1/6.0 + a2*a2/10.0)*Lap_K(r,s,i,0,2)(0,0))),
                  S2*(-c*(K12(0,0) + (a1*a1/6.0 + a2*a2/10.0)*Lap_K(r,s,i,1,2)(0,0)));
    } else {
        result.setZero();
    }
    return result;
}

// Add new tensor operation to calculate D_K
Matrix3d D_K(const Vector3d& r, double s, int l, int i, int j, int k) {
    return 0.5*(DD_J(r,s,l,k,i,j) + DD_J(r,s,l,j,i,k));
}

// Helper for M23 calculation - move this BEFORE con_M23_row
double sum_levi_DK(const Vector3d& r, double s, int i, int j1, int j2, double a2) {
    double sum = 0;
    for(int l = 0; l < 3; l++) {
        for(int m = 0; m < 3; m++) {
            if(l != m && i != l && i != m) {
                Matrix3d dk = D_K(r,s,l,m,j1,j2);
                Matrix3d dlk = DLap_K(r,s,l,m,j1,j2);
                sum += levi(i,l,m) * (dk(0,0) + (a2*a2/6.0)*dlk(0,0));
            }
        }
    }
    return sum;
}

// Now con_M23_row can use sum_levi_DK
Vector5d con_M23_row(const Vector3d& r, double s, double a1, double a2, int i, double c, double mu) {
    Vector5d result;
    if(s > 1e-10) {
        Matrix3d sum = Matrix3d::Zero();
        for(int l = 0; l < 3; l++) {
            for(int m = 0; m < 3; m++) {
                if(l != m && i != l && i != m) {
                    sum += levi(i,l,m) * (D_K(r,s,l,m,0,0) + (a2*a2/6.0)*DLap_K(r,s,l,m,0,0));
                }
            }
        }
        double A = -0.5*c*sum(0,0);
        
        sum = Matrix3d::Zero();
        for(int l = 0; l < 3; l++) {
            for(int m = 0; m < 3; m++) {
                if(l != m && i != l && i != m) {
                    sum += levi(i,l,m) * (D_K(r,s,l,m,1,1) + (a2*a2/6.0)*DLap_K(r,s,l,m,1,1));
                }
            }
        }
        double B = -0.5*c*sum(0,0);

        // Calculate other components
        double C01 = -0.5*c*sum_levi_DK(r,s,i,0,1,a2);
        double C02 = -0.5*c*sum_levi_DK(r,s,i,0,2,a2);
        double C12 = -0.5*c*sum_levi_DK(r,s,i,1,2,a2);
        
        result << (HS3P1*A + HS3M1*B),
                  S2*C01,
                  (HS3M1*A + HS3P1*B),
                  S2*C02,
                  S2*C12;
    } else {
        result.setZero();
    }
    return result;
}

// Fix D_R to handle scalar operations properly 
Matrix3d D_R(const Vector3d& r, double s, int l, int i, int j) {
    Matrix3d result = Matrix3d::Zero();
    for(int k = 0; k < 3; k++) {
        for(int m = 0; m < 3; m++) {
            if(k != m && j != k && j != m) {
                Matrix3d DJ = DD_J(r,s,l,k,i,m);
                double lev = static_cast<double>(levi(j,k,m));
                result += lev * DJ;
            }
        }
    }
    return -0.5 * result;
}

// Fix M22 block calculation
Matrix3d M22_calc(const Vector3d& r, double s, double a1, double a2, double c, double mu) {
    Matrix3d result = Matrix3d::Zero();
    if(s > 1e-10) {
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                double sum = 0.0;
                for(int k = 0; k < 3; k++) {
                    for(int l = 0; l < 3; l++) {
                        if(k != l && i != k && i != l) {
                            Matrix3d DR = D_R(r,s,k,l,j);
                            sum += levi(i,k,l) * DR(0,0);
                        }
                    }
                }
                result(i,j) = 0.5 * c * sum;
            }
        }
    } else {
        result = Matrix3d::Identity() / (8.0 * PI * mu * pow(a1,3));
    }
    return result;
}

// Fix M33 block calculation
void calc_M33_block(Matrix<double,5,5>& M33_temp, const Vector3d& r, double s, 
                    double a1, double a2, double c, double mu) {
    for(int p = 0; p < 5; p++) {
        for(int q = 0; q < 5; q++) {
            int i1 = COND_IDX[p][0], i2 = COND_IDX[p][1];
            int j1 = COND_IDX[q][0], j2 = COND_IDX[q][1];
            
            Matrix3d dk1 = D_K(r,s,j1,i1,i2,j2);
            Matrix3d dk2 = D_K(r,s,i1,j1,i2,j2);
            Matrix3d dlk1 = DLap_K(r,s,j1,i1,i2,j2);
            Matrix3d dlk2 = DLap_K(r,s,i1,j1,i2,j2);
            
            M33_temp(p,q) = -0.5*c*((dk1(0,0) + dk2(0,0)) + 
                                   (a1*a1 + a2*a2)/10.0 * (dlk1(0,0) + dlk2(0,0)));
        }
    }
}

// Main Minfinity computation with OpenMP optimization
MatrixXd computeMinfinity(const std::vector<Vector3d>& positions, const std::vector<double>& radii, double mu) {
    if (positions.size() != radii.size()) {
        throw std::invalid_argument("Number of positions must match number of radii");
    }
    if (mu <= 0) {
        throw std::invalid_argument("Viscosity must be positive");
    }
    
    int n = positions.size();
    int matrix_size = 11*n; // Full size including rotation and strain blocks
    MatrixXd M = MatrixXd::Zero(matrix_size, matrix_size);
    
    double c = 1.0/(8.0 * PI * mu);
    
    // Pre-allocate thread-local matrices to avoid repeated allocations
    std::vector<MatrixXd> thread_local_Ms;
    int num_threads = 1;
    
    #ifdef _OPENMP
    num_threads = omp_get_max_threads();
    #endif
    
    thread_local_Ms.resize(num_threads, MatrixXd::Zero(matrix_size, matrix_size));
    
    #pragma omp parallel
    {
        int thread_id = 0;
        #ifdef _OPENMP
        thread_id = omp_get_thread_num();
        #endif
        
        MatrixXd& local_M = thread_local_Ms[thread_id];
        
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n; i++) {
            if (radii[i] <= 0) {
                #pragma omp critical
                {
                    throw std::invalid_argument("Particle radii must be positive");
                }
            }
            
            for (int j = 0; j < n; j++) {
                // Get block indices for all submatrices
                int A_i = 3*i, A_j = 3*j;           // M11 block
                int Bt_i = 3*i, Bt_j = 3*n + 3*j;   // M12 block
                int C_i = 3*n + 3*i, C_j = 3*n + 3*j; // M22 block
                int G_i = 3*i, G_j = 6*n + 5*j;     // M13 block
                int H_i = 3*n + 3*i, H_j = 6*n + 5*j; // M23 block
                int M_i = 6*n + 5*i, M_j = 6*n + 5*j; // M33 block
                
                if (i == j) {
                    // Diagonal blocks
                    local_M.block<3,3>(A_i, A_j) = Matrix3d::Identity() / (6.0 * PI * mu * radii[i]);
                    local_M.block<3,3>(C_i, C_j) = Matrix3d::Identity() / (8.0 * PI * mu * pow(radii[i],3));
                    // Diagonal M33 block
                    for(int k = 0; k < 5; k++) {
                        local_M(M_i+k, M_j+k) = 1.0/(20.0/3.0 * PI * mu * pow(radii[i],3));
                    }
                } else {
                    Vector3d r = positions[i] - positions[j];
                    double s = r.norm();
                    
                    // Handle close particles
                    if (s > 1e-8 && 2*s/(radii[i]+radii[j]) < 2.001) {
                        double ss_out = 2.001 * (radii[i]+radii[j]) / 2.0;
                        r *= ss_out/s;
                        s = ss_out;
                    }
                    
                    if (s > 1e-8) {
                        double a1 = radii[i], a2 = radii[j];
                        
                        // M11 block
                        local_M.block<3,3>(A_i, A_j) = c * (J_tensor(r,s) + 
                            (a1*a1 + a2*a2)/6.0 * Lap_J(r,s));
                        
                        // M12 block
                        local_M.block<3,3>(Bt_i, Bt_j) = c * (R_tensor(r,s) + 
                            a1*a1/6.0 * Lap_R(r,s));
                            
                        // M22 block - Fixed calculation
                        local_M.block<3,3>(C_i, C_j) = M22_calc(r, s, a1, a2, c, mu);
                        
                        // M13 block (G)
                        for(int k = 0; k < 3; k++) {
                            local_M.block<1,5>(G_i+k, G_j) = con_M13_row(r, s, a1, a2, k, c, mu).transpose();
                        }
                        
                        // M23 block (H)
                        for(int k = 0; k < 3; k++) {
                            local_M.block<1,5>(H_i+k, H_j) = con_M23_row(r, s, a1, a2, k, c, mu).transpose();
                        }
                        
                        // M33 block
                        Matrix<double,5,5> M33_temp = Matrix<double,5,5>::Zero();
                        calc_M33_block(M33_temp, r, s, a1, a2, c, mu);
                        local_M.block<5,5>(M_i, M_j) = COND_E * M33_temp * COND_E;
                    } else {
                        // Diagonal case for M33
                        for(int k = 0; k < 5; k++) {
                            local_M(M_i+k, M_j+k) = 1.0/(20.0/3.0 * PI * mu * pow(radii[i],3));
                        }
                    }
                }
            }
        }
    }
    
    // Combine thread-local results
    for (const auto& local_M : thread_local_Ms) {
        M += local_M;
    }
    
    return 0.5 * (M + M.transpose());
}

// Compute Rinfinity - the inverse of Minfinity
MatrixXd computeRinfinity(const std::vector<Vector3d>& positions, const std::vector<double>& radii, double mu) {
    // First compute the mobility matrix
    MatrixXd M = computeMinfinity(positions, radii, mu);
    
    // Make it symmetric to ensure numerical stability
    M = 0.5 * (M + M.transpose());
    
    // Use Eigen's direct solver for inversion
    // For larger matrices, we might want to use LU decomposition for better performance
    int n = positions.size();
    if (n > 100) { // Threshold for using LU decomposition
        return M.partialPivLu().inverse();
    } else {
        return M.inverse(); // For smaller matrices, direct inverse is fine
    }
}

#ifdef _OPENMP
// OpenMP thread control functions
int get_max_threads() {
    return omp_get_max_threads();
}

void set_num_threads(int num_threads) {
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
}

int get_num_threads() {
    int num_threads = 0;
    #pragma omp parallel
    {
        #pragma omp single
        num_threads = omp_get_num_threads();
    }
    return num_threads;
}
#else
// Fallback implementations when OpenMP is not available
int get_max_threads() {
    return 1;  // Only one thread available without OpenMP
}

void set_num_threads(int num_threads) {
    // Do nothing - OpenMP not available
}

int get_num_threads() {
    return 1;  // Only one thread available without OpenMP
}
#endif
