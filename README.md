# StokesKit

<p align="center">
    <br/>
    <img src="images/logo.png" alt="StokesKit Logo" width="300"/>
    <br/>
</p>


StokesKit is a high-performance C++/Python library for computing hydrodynamic interactions between particles in Stokes flow. It provides efficient implementations of scalar resistance functions, mobility matrices, and exact two-body resistance matrices for low Reynolds number hydrodynamics.

> **⚠️ WARNING: Development Status**  
> StokesKit is currently in active development. The API may change, and some implementations may not be fully correct or optimized. Use with caution in production environments. Feedback and bug reports are welcome!

## Features

- **Scalar Resistance Functions**: XA, YA, YB, XC, YC functions based on Jeffrey & Onishi (1984)
- **Mobility Matrices**: Fast computation of many-body mobility (Minfinity) matrices
- **Resistance Matrices**: Exact two-body resistance matrices (R2Bexact)
  - Note: The rate-of-strain components in the R2Bexact matrix are not yet fully implemented but will be available soon
- **Optimized Performance**: Intelligent caching system for function evaluations
- **Parallelization**: OpenMP support for multi-threading
- **Memory Profiling**: Optional memory usage tracking

## Installation

### Prerequisites

- C++ compiler with C++17 support
- Python 3.6 or higher
- NumPy
- OpenMP

### Installing from source

```bash
# Clone the repository
git clone https://github.com/dimos99/stokeskit.git
cd stokeskit

# Install with pip (standard)
pip install .

# Install with profiling enabled
pip install . --install-option="--enable-profiling"

# Install in development mode
pip install -e .

# Install with debug symbols
pip install . --install-option="--debug"
```

## Usage

### Python Interface

```python
import numpy as np
import stokeskit as sk

# Computing scalar resistance functions
s_values = np.linspace(2.001, 10.0, 100)  # dimensionless separations
l_values = np.ones_like(s_values)         # size ratios (1 = equal sized)

# Compute all resistance functions at once
result = sk.compute_resistance_functions(s_values, l_values)
xa11 = result["XA11"]
xa12 = result["XA12"]

# Computing mobility matrix for a collection of particles
positions = np.array([
    [0, 0, 0],
    [0, 0, 3],
    [0, 3, 0]
])
radii = np.array([1.0, 1.0, 1.0])
viscosity = 1.0

# Compute mobility matrix
M = sk.Minfinity.computeMinfinityArray(positions, radii, viscosity)

# Compute resistance matrix
R = sk.Minfinity.computeRinfinityArray(positions, radii, viscosity)

# Compute exact two-body resistance matrix
R2B = sk.R2Bexact.generate(positions, radii, viscosity)
```

### OpenMP Parallelization

```python
# Check if OpenMP is available
if sk.Minfinity.has_openmp():
    # Get and set the number of threads
    max_threads = sk.Minfinity.get_max_threads()
    sk.Minfinity.set_num_threads(4)  # Set to 4 threads
```

### Individual Resistance Functions

```python
# Computing individual resistance functions
# Dimensionless separation and radius ratio
s = 3.0  # s = r/(a1+a2)
l = 1.0  # l = a2/a1

# Compute scalar resistance functions
xa11 = sk.XA.XA11(s, l)
ya11 = sk.YA.YA11(s, l)
yb11 = sk.YB.YB11(s, l)
xc11 = sk.XC.XC11(s, l)
yc11 = sk.YC.YC11(s, l)

# Clear caches to free memory
sk.XA.clear_cache()
sk.YA.clear_cache()
sk.YB.clear_cache()
sk.XC.clear_cache()
sk.YC.clear_cache()
```

## Examples

### 1. Computing Resistance Functions for a Range of Separations
```python
import numpy as np
import stokeskit as sk

# Configure calculation parameters
lubr_cutoff = 2.01  # Lubrication regime cutoff
nf_cutoff = 4.0     # Near-field to far-field transition
maxIter = 300       # Maximum iterations for convergence
rtol = 1e-5         # Relative tolerance
atol = 1e-8         # Absolute tolerance

# Create a log-spaced range of separation parameters
ksis = np.logspace(-5, 2, 1000)  # ξ = s-2 (dimensionless surface separation)
l_value = 0.5                    # λ = a₂/a₁ (size ratio)
ss = 2 + ksis                    # s = r/(a₁+a₂) (dimensionless center-to-center distance)

# Compute all resistance functions at once (vectorized)
resistance = sk.compute_resistance_functions(
    ss, l_value * np.ones_like(ss),
    lubr_cutoff=lubr_cutoff,
    cutoff=nf_cutoff,
    maxIter=maxIter,
    rtol=rtol,
    atol=atol
)

# Extract specific functions
xa11 = resistance["XA11"]
xa12 = resistance["XA12"]
ya11 = resistance["YA11"]
ya12 = resistance["YA12"]
yb11 = resistance["YB11"]
yb12 = resistance["YB12"]
xc11 = resistance["XC11"]
xc12 = resistance["XC12"]
yc11 = resistance["YC11"]
yc12 = resistance["YC12"]

# Example of analyzing the lubrication behavior
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.loglog(ksis, xa11, label="XA11")
plt.loglog(ksis, xa12, label="XA12")
plt.xlabel(r"$\xi = s-2$")
plt.ylabel("Resistance Functions")
plt.grid(True, alpha=0.3)
plt.legend()
plt.title(f"Resistance Functions (λ={l_value})")
plt.savefig("resistance_functions.png")
```

### 2. Computing the Mobility Matrix with Minfinity

```python
import numpy as np
import stokeskit as sk
import time

# Set up a system of particles
mu = 1.0  # Fluid viscosity

# Create a configuration of particles in a triangle
positions = np.array([
    [0.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
    [1.5, 2.6, 0.0]
])
radii = np.array([1.0, 1.0, 1.0])

# Set OpenMP threads if available
if hasattr(sk.Minfinity, 'has_openmp') and sk.Minfinity.has_openmp():
    max_threads = sk.Minfinity.get_max_threads()
    print(f"Using {max_threads} OpenMP threads")
    sk.Minfinity.set_num_threads(max_threads)

# Compute the mobility matrix
start = time.time()
M = sk.Minfinity.computeMinfinityArray(positions, radii, mu)
elapsed = time.time() - start
print(f"Computed mobility matrix in {elapsed:.4f} seconds")

# Extract the translation-translation block (first 9x9 for 3 particles)
M_tt = M[:9, :9]
print("Translation-Translation block of mobility matrix:")
print(M_tt)

# Compute the resistance matrix directly (instead of inverting M)
R = sk.Minfinity.computeRinfinityArray(positions, radii, mu)
print("Resistance matrix computed directly:")
print(R[:9, :9])  # Translation-Translation block
```

### 3. Generating Exact Two-Body Resistance Matrix with R2Bexact

```python
import numpy as np
import stokeskit as sk

# Set up two particles along the z-axis
positions = np.array([
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 3.0]  # Separation of 3 units along z-axis
])
radii = np.array([1.0, 1.0])  # Equal-sized particles
mu = 1.0  # Fluid viscosity

# Generate the exact two-body resistance matrix
R = sk.R2Bexact.generate(positions, radii, mu)

# The full matrix has size 22x22 for two particles
# Each particle contributes 11 degrees of freedom:
# - 3 for translation
# - 3 for rotation
# - 5 for rate of strain (coming soon)

# Extract the translation-translation block (6x6 for 2 particles)
R_tt = R[:6, :6]
print("Translation-Translation block of R2Bexact:")
print(R_tt)

# For unequal particles, try different size ratios
size_ratios = [0.5, 1.0, 2.0]
for ratio in size_ratios:
    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 3.0 * (1+ratio)]  # Adjust separation based on size
    ])
    radii = np.array([1.0, ratio])
    
    # Generate matrix and extract diagonal elements
    R = sk.R2Bexact.generate(positions, radii, mu)
    diag = np.diag(R[:6, :6])
    print(f"Size ratio {ratio}, diagonal elements: {diag}")
```

<!-- ### 4. Validation and Benchmarking (not ready yet)

StokesKit includes a comprehensive validation script (`check_code.py`) that performs:

1. Comparison of resistance functions with analytical expressions
2. Performance benchmarking of vectorized vs. individual calculations
3. Testing Minfinity with varying numbers of particles and OpenMP threads
4. Validation of R2Bexact against far-field (Minfinity) approximations
5. Generation of visualization plots for resistance functions

```bash
# Run the full validation and benchmarking suite
python check_code.py

# This will generate several plots:
# - comparison_l*.png - Comparing vectorized vs individual calculations
# - XA_YA_functions.png - XA and YA functions for different λ values
# - YB_functions.png - YB functions for different λ values
# - XC_YC_functions.png - XC and YC functions for different λ values
# - minfinity_performance*.png - Minfinity scaling performance
# - r2bexact_*.png - R2Bexact validation and performance
``` -->

## Performance Considerations

### Vectorized vs. Individual Calculations

StokesKit offers two methods to calculate resistance functions:

1. **Individual functions** (`XA.XA11(s, l)`, etc.): Best for single-point evaluations
2. **Vectorized calculations** (`compute_resistance_functions(ss, ls)`): Significantly faster for batch calculations

The vectorized implementation uses OpenMP parallelization and optimized algorithms to compute all resistance functions simultaneously for arrays of inputs.

### First-time Calculation Overhead

When scalar resistance functions are calculated for the first time, especially for complex parameter ranges, the computation may take significant time due to the series expansions. StokesKit implements two separate caching mechanisms to mitigate this overhead:

1. **Persistent Cache (.cache folder)**: StokesKit creates a compressed cache inside a `.cache` folder for difficult calculations. This cache persists between program executions, allowing subsequent runs to reuse previously computed values.

2. **In-memory Cache**: During program execution, StokesKit maintains an on-the-fly caching mechanism that stores computed values in memory. This cache is faster than the persistent cache but doesn't persist after the program ends.

Both caching mechanisms dramatically improve performance for repeated calculations:

```python
# Example demonstrating cache speedup
import time
import stokeskit as sk

# First calculation (might be slow - both caches are empty)
start = time.time()
result1 = sk.compute_resistance_functions(s_values, l_values)
print(f"First calculation: {time.time() - start:.2f} seconds")

# Second calculation with same parameters (much faster - using in-memory cache)
start = time.time()
result2 = sk.compute_resistance_functions(s_values, l_values)
print(f"Second calculation: {time.time() - start:.2f} seconds")

# If you restart the program, the persistent .cache will be used
# so it will still be faster than the initial calculation
```

### Caching System

The in-memory caching system can be managed to free up memory when no longer needed:

```python
# Clear all in-memory caches
sk.XA.clear_cache()
sk.YA.clear_cache()
sk.YB.clear_cache()
sk.XC.clear_cache()
sk.YC.clear_cache()
```

### Memory Profiling

When built with memory profiling enabled, you can track memory usage:

```python
# Print memory usage statistics (if compiled with ENABLE_PROFILING)
sk.print_memory_usage()
sk.print_profile_stats()
```

### OpenMP Scaling

The Minfinity and vectorized resistance calculations scale efficiently with additional OpenMP threads. For optimal performance:

1. Use vectorized calculations for large parameter sweeps
2. Set an appropriate thread count for your hardware

## Theoretical Background

The hydrodynamic interactions are based on the resistance formulation described in:

1. Jeffrey, D.J. and Onishi, Y. (1984), "Calculation of the resistance and mobility functions for two unequal rigid spheres in low-Reynolds-number flow", Journal of Fluid Mechanics, 139, 261–290.

2. Kim, S. and Karrila, S.J. (1991), "Microhydrodynamics: Principles and Selected Applications", Butterworth-Heinemann.

3. Happel, J. and Brenner, H. (1983), "Low Reynolds Number Hydrodynamics: With Special Applications to Particulate Media", Springer Science & Business Media.

4. Brady, J.F. and Bossis, G. (1988), "Stokesian Dynamics", Annual Review of Fluid Mechanics, 20, 111–157.


StokesKit implements three regimes for resistance functions:
- **Lubrication regime** (s < lubr_cutoff): Asymptotic expressions for very close particles
- **Near-field regime** (lubr_cutoff ≤ s < cutoff): Series expansion for moderate separations
- **Far-field regime** (s ≥ cutoff): Series expansion for well-separated particles

## License

[![License: EUPL](https://img.shields.io/badge/License-EUPL%201.2-blue.svg)](https://joinup.ec.europa.eu/collection/eupl/eupl-text-11-12)

Licensed under the EUPL.

For the full license text, please see: https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12

## Contributing

Contributions to StokesKit are welcome! Please feel free to submit a Pull Request.

## Citation

If you use StokesKit in your research, please cite it as:

```bibtex
@software{aslanis_stokeskit_2025,
  author       = {Aslanis, Dimos},
  title        = {StokesKit: A high-performance library for hydrodynamic interactions in Stokes flow},
  year         = {2025},
  url          = {https://github.com/dimos99/stokeskit},
  affiliation  = {Utrecht University}
}
```

## Acknowledgements

This work is funded from the European Union's Horizon Europe Framework Programme (HORIZON) under the Marie Skłodowska-Curie Grant Agreement (GA) Nº: 101120301
For more information, please visit the [project website](https://cocogel.iesl.forth.gr).

The project was done in Utrecht University, under the supervision of Prof. Joost de Graaf.

The logo was designed by [kah3nee](https://github.com/kah3nee), inspired by [The Strokes](https://www.thestrokes.com/) band logo.

## Contact

For questions or support, please contact the author at: `d.aslanis@uu.nl`