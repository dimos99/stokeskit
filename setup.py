from setuptools import setup, Extension
import pybind11
import sys
import platform
import numpy as np
from setuptools.command.build_ext import build_ext
import os

# Add command line options
enable_profiling = '--enable-profiling' in sys.argv
enable_debug = '--debug' in sys.argv

if enable_profiling:
    sys.argv.remove('--enable-profiling')
if enable_debug:
    sys.argv.remove('--debug')

# Define compiler flags and preprocessor definitions
extra_compile_args = ['-std=c++17', '-fvisibility=default']
extra_link_args = []  # Initialize empty, will be set per platform

# Add OpenMP flags based on platform
if sys.platform == 'darwin':  # macOS
    extra_compile_args.extend(['-Xpreprocessor', '-fopenmp'])
    extra_link_args = ['-undefined', 'dynamic_lookup']  # macOS specific
    if os.path.exists('/opt/homebrew/opt/libomp/lib/'):  # Apple Silicon
        extra_link_args.extend(['-L/opt/homebrew/opt/libomp/lib/', '-lomp'])
        extra_compile_args.extend(['-I/opt/homebrew/opt/libomp/include/'])
    elif os.path.exists('/usr/local/opt/libomp/lib/'):  # Intel Mac
        extra_link_args.extend(['-L/usr/local/opt/libomp/lib/', '-lomp'])
        extra_compile_args.extend(['-I/usr/local/opt/libomp/include/'])
    else:
        raise Exception("OpenMP not found! Please install it with 'brew install libomp'")
elif sys.platform == 'linux':  # Linux
    extra_compile_args.append('-fopenmp')
    extra_link_args.extend(['-fopenmp', '-lz'])  # Add zlib
else:  # Windows
    extra_compile_args.append('/openmp')

define_macros = []

if enable_profiling:
    define_macros.append(('ENABLE_PROFILING', '1'))

if enable_debug:
    define_macros.append(('DEBUG', '1'))
    if sys.platform == 'darwin' or sys.platform == 'linux':
        extra_compile_args.extend(['-g', '-O0', '-Wall', '-Wextra'])
    else:  # Windows
        extra_compile_args.extend(['/Od', '/Zi', '/Wall'])
else:
    # Enhanced Release mode flags
    if sys.platform == 'darwin' or sys.platform == 'linux':
        extra_compile_args.extend([
            '-O3',                    # Maximum optimization
            '-march=native',          # CPU-specific optimizations
            '-ffast-math',            # Aggressive floating-point optimizations
            '-flto',                  # Link-time optimization
            '-fuse-linker-plugin',    # Better LTO
            '-funroll-loops',         # Loop unrolling
            '-mavx2',                 # Enable AVX2 instructions if available
            '-mfma',                  # Enable FMA instructions if available
            '-DNDEBUG'
        ])
        extra_link_args.extend([
            '-flto',                  # Link-time optimization
            '-fuse-linker-plugin'     # Better LTO
        ])
    else:  # Windows
        extra_compile_args.extend([
            '/O2',                    # Maximum optimization
            '/GL',                    # Whole program optimization
            '/fp:fast',              # Fast floating-point model
            '/arch:AVX2',            # Enable AVX2 instructions
            '/DNDEBUG'
        ])
        extra_link_args = ['/LTCG']  # Link-time code generation

# Create the extension module
ext_modules = [
    Extension(
        "stokeskit",
        ["src/stokeskit.cpp",
         "src/common.cpp",
         "src/XA.cpp",
         "src/YA.cpp",
         "src/YB.cpp",
         "src/XC.cpp",
         "src/YC.cpp"],         
        include_dirs=[pybind11.get_include(), np.get_include(), 'src'],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        define_macros=define_macros,
        language='c++'
    )
]

setup(
    name="stokeskit",
    version="0.1.0",
    author="Dimos Aslanis",
    author_email="d.aslanis@uu.nl",
    description="Resistance scalar functions for two-body hydrodynamic interactions",
    long_description="",
    ext_modules=ext_modules,
    setup_requires=["pybind11>=2.6.0"],
    install_requires=["pybind11>=2.6.0"],
    python_requires=">=3.6",
)
