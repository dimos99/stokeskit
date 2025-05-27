from setuptools import setup, Extension
import pybind11
import sys
import platform
import numpy as np
from setuptools.command.build_ext import build_ext
import os
import urllib.request
import zipfile
import shutil

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

# Initialize include_dirs early
include_dirs = []

# Detect architecture
is_arm = platform.machine() in ('arm64', 'aarch64')

# Function to check if a file exists in any of the directories
def find_file(filename, search_paths):
    for path in search_paths:
        file_path = os.path.join(path, filename)
        if os.path.exists(file_path):
            return path
    return None

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
    
    # Check for zlib headers
    zlib_search_paths = [
        '/usr/include',
        '/usr/local/include',
        '/opt/local/include',
        '/home/dimos/.conda/envs/jaxsd/include',
    ]
    
    zlib_include_path = find_file('zlib.h', zlib_search_paths)
    if zlib_include_path:
        print(f"Found zlib.h in {zlib_include_path}")
        if zlib_include_path not in include_dirs:
            include_dirs.append(zlib_include_path)
    else:
        print("ERROR: zlib.h not found! Please install zlib development headers.")
        print("On Ubuntu/Debian: sudo apt-get install zlib1g-dev")
        print("On CentOS/RHEL: sudo yum install zlib-devel")
        print("On Conda: conda install -c anaconda zlib")
        sys.exit(1)
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
        # Common optimization flags for all architectures
        arch_specific_flags = [
            '-O3',                    # Maximum optimization
            '-march=native',          # CPU-specific optimizations
            '-ffast-math',            # Aggressive floating-point optimizations
            '-funroll-loops',         # Loop unrolling
            '-DNDEBUG'
        ]
        
        # LTO flags - don't use fuse-linker-plugin on ARM
        if not is_arm:
            arch_specific_flags.append('-flto')
            extra_link_args.append('-flto')
            
            # x86 specific flags
            arch_specific_flags.extend([
                '-mavx2',                 # Enable AVX2 instructions (x86 only)
                '-mfma'                   # Enable FMA instructions (x86 only)
            ])
            
            # Only add fuse-linker-plugin on x86
            arch_specific_flags.append('-fuse-linker-plugin')
            extra_link_args.append('-fuse-linker-plugin')
        else:
            # ARM specific LTO flag without plugin
            arch_specific_flags.append('-flto')
            extra_link_args.append('-flto')
            
        extra_compile_args.extend(arch_specific_flags)
    else:  # Windows
        extra_compile_args.extend([
            '/O2',                    # Maximum optimization
            '/GL',                    # Whole program optimization
            '/fp:fast',              # Fast floating-point model
            '/arch:AVX2',            # Enable AVX2 instructions
            '/DNDEBUG'
        ])
        extra_link_args = ['/LTCG']  # Link-time code generation

# Ensure third_party directory exists
third_party_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'third_party')
if not os.path.exists(third_party_dir):
    os.makedirs(third_party_dir)

# Download and extract Eigen if not present - fixes to ensure proper structure
eigen_dir = os.path.join(third_party_dir, 'eigen')
eigen_include_dir = os.path.join(eigen_dir, 'include')

if not os.path.exists(eigen_include_dir):
    print("Downloading Eigen library...")
    eigen_url = "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip"
    zip_path = os.path.join(third_party_dir, "eigen.zip")
    
    # Create directories
    os.makedirs(eigen_dir, exist_ok=True)
    os.makedirs(eigen_include_dir, exist_ok=True)
    
    # Download the zip file
    urllib.request.urlretrieve(eigen_url, zip_path)
    
    # Extract the zip file
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(third_party_dir)
    
    # Get the extracted directory path
    extracted_dir = os.path.join(third_party_dir, "eigen-3.4.0")
    
    # Create Eigen directory structure
    if os.path.exists(extracted_dir):
        # Move Eigen header files to include/Eigen
        eigen_headers_dir = os.path.join(eigen_include_dir, "Eigen")
        os.makedirs(eigen_headers_dir, exist_ok=True)
        
        # Copy all Eigen headers (recursively) to include/Eigen
        for item in os.listdir(os.path.join(extracted_dir, "Eigen")):
            src = os.path.join(extracted_dir, "Eigen", item)
            dst = os.path.join(eigen_headers_dir, item)
            if os.path.isdir(src):
                shutil.copytree(src, dst)
            else:
                shutil.copy2(src, dst)
        
        # Copy unsupported directory if it exists
        unsupported_src = os.path.join(extracted_dir, "unsupported")
        if os.path.exists(unsupported_src):
            unsupported_dst = os.path.join(eigen_include_dir, "unsupported")
            shutil.copytree(unsupported_src, unsupported_dst)
        
        # Clean up the extracted directory
        shutil.rmtree(extracted_dir)
    
    # Clean up zip file
    os.remove(zip_path)
    print("Eigen library installed successfully.")

# Include directories - add Eigen include path and potential zlib paths
include_dirs += [
    pybind11.get_include(), 
    np.get_include(), 
    'src',
    eigen_include_dir,  # Changed to point directly to the include directory
]

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
         "src/YC.cpp",
         "src/Minfinity.cpp",
         "src/R2Bexact.cpp"],
        include_dirs=include_dirs,
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
