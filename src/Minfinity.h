#ifndef MINFINITY_H
#define MINFINITY_H

#include <vector>
#include <array>
#include <Eigen/Dense>

// Type definition for 5D vector
typedef Eigen::Matrix<double, 5, 1> Vector5d;

// Main function to compute the Minfinity matrix
Eigen::MatrixXd computeMinfinity(
    const std::vector<Eigen::Vector3d>& positions, 
    const std::vector<double>& radii, 
    double mu);

// Function to compute Rinfinity (inverse of Minfinity)
Eigen::MatrixXd computeRinfinity(
    const std::vector<Eigen::Vector3d>& positions, 
    const std::vector<double>& radii, 
    double mu);

// Helper functions for tensor operations - exposed for testing purposes
int levi(int i, int j, int k);
Eigen::Matrix3d J_tensor(const Eigen::Vector3d& r, double s);
Eigen::Matrix3d Lap_J(const Eigen::Vector3d& r, double s);
Eigen::Matrix3d D_J(const Eigen::Vector3d& r, double s, int l, int i, int j);
Eigen::Matrix3d R_tensor(const Eigen::Vector3d& r, double s);
Eigen::Matrix3d K_tensor(const Eigen::Vector3d& r, double s, int i, int j, int k);
Eigen::Matrix3d DD_J(const Eigen::Vector3d& r, double s, int m, int l, int i, int j);
Eigen::Matrix3d DLap_J(const Eigen::Vector3d& r, double s, int k, int i, int j);
Eigen::Matrix3d Lap_R(const Eigen::Vector3d& r, double s);
Eigen::Matrix3d Lap_K(const Eigen::Vector3d& r, double s, int i, int j, int k);
Eigen::Matrix3d D_K(const Eigen::Vector3d& r, double s, int l, int i, int j, int k);
Eigen::Matrix3d DLap_K(const Eigen::Vector3d& r, double s, int l, int i, int j, int k);
Eigen::Matrix3d D_R(const Eigen::Vector3d& r, double s, int l, int i, int j);

// Matrix contraction helpers
Vector5d con_M13_row(const Eigen::Vector3d& r, double s, double a1, double a2, int i, double c, double mu);
double sum_levi_DK(const Eigen::Vector3d& r, double s, int i, int j1, int j2, double a2);
Vector5d con_M23_row(const Eigen::Vector3d& r, double s, double a1, double a2, int i, double c, double mu);
void calc_M33_block(Eigen::Matrix<double,5,5>& M33_temp, const Eigen::Vector3d& r, double s, 
                   double a1, double a2, double c, double mu);

#endif // MINFINITY_H
