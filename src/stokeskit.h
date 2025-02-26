#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

// Import all submodule headers
#include "XA.h"
#include "YA.h"
#include "YB.h"
#include "XC.h"
#include "YC.h"
// Future headers will go here (YA, XB, YB, XC, YC, etc.)

namespace py = pybind11;

// Forward declarations of module creation functions
void init_XA(py::module& m);
void init_YA(py::module& m);
void init_YB(py::module& m);
void init_XC(py::module& m);
void init_YC(py::module& m);
// Future module init declarations will go here
