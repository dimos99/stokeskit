#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>  // Add for Eigen support

// Import all submodule headers
#include "XA.h"
#include "YA.h"
#include "YB.h"
#include "XC.h"
#include "YC.h"
#include "Minfinity.h"  // Add Minfinity header

namespace py = pybind11;

// Forward declarations of module creation functions
void init_XA(py::module& m);
void init_YA(py::module& m);
void init_YB(py::module& m);
void init_XC(py::module& m);
void init_YC(py::module& m);
void init_Minfinity(py::module& m);  // Add Minfinity initialization function
