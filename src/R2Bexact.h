#pragma once

#include <Eigen/Dense>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include "common.h"
#include "XA.h"
#include "YA.h"
#include "YB.h"
#include "XC.h"
#include "YC.h"

using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::MatrixXd;
using Eigen::VectorXd;

const double PI = 3.141592653589793;
const double S2 = std::sqrt(2.0);
const double S3 = std::sqrt(3.0);
const double HS3P1 = 0.5*(S3+1);
const double HS3M1 = 0.5*(S3-1);

// Condensed index definitions
extern const std::array<std::array<int,2>,5> COND_IDX;
extern const Eigen::MatrixXd COND_E;

// Levi-Civita tensor is already declared in Minfinity.h
// Just declare it here without the "extern" keyword so it's treated as a forward declaration
int levi(int i, int j, int k);

// Unit displacement tensors
double L1(const Vector3d& d, int i, int j);
double L2(const Vector3d& d, int i, int j);
double L3(const Vector3d& d, int i, int j);
double L4(const Vector3d& d, int i, int j, int k);
double L5(const Vector3d& d, int i, int j, int k);
double L6(const Vector3d& d, int i, int j, int k);
double L7(const Vector3d& d, int i, int j, int k, int l);
double L8(const Vector3d& d, int i, int j, int k, int l);
double L9(const Vector3d& d, int i, int j, int k, int l);

// Helper functions for the R2Bexact submatrix elements
double Af(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i, int j);
double Bf(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i, int j);
double Cf(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i, int j);
double Gf(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i, int j, int k);
double Hf(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i, int j, int k);
double Mf(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i, int j, int k, int l);

// Condensed submatrix helper functions
Eigen::VectorXd con_Gf_row(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i);
Eigen::VectorXd con_Hf_row(int gamma, const Eigen::Vector3d& d, double lambda, double s, int i);

// The main function to generate the R2Bexact resistance matrix
Eigen::MatrixXd generateR2Bexact(
    const std::vector<Eigen::Vector3d>& positions,
    const std::vector<double>& radii,
    double mu = 1.0,
    double cutoff_factor = 2.0
);

// Helper function for close particle detection
std::tuple<std::vector<std::pair<int, int>>, std::vector<Eigen::Vector3d>, 
          std::vector<double>, std::vector<double>> 
findCloseParticles(
    const std::vector<Eigen::Vector3d>& positions,
    const std::vector<double>& radii,
    double cutoff_factor
);

// Submatrix index helpers
struct SubmatrixIndices {
    int A_start_row, A_start_col;
    int B_start_row, B_start_col;
    int B21_start_row, B21_start_col;
    int C_start_row, C_start_col;
    int G_start_row, G_start_col;
    int G21_start_row, G21_start_col;
    int H_start_row, H_start_col;
    int H21_start_row, H21_start_col;
    int M_start_row, M_start_col;
};

SubmatrixIndices calculateSubmatrixIndices(int a1_index, int a2_index, int num_spheres);
