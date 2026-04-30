# How ecbuild Finds the Intel MKL Library

This document explains the mechanism used by `ecbuild/cmake/FindMKL.cmake` to locate the Intel Math Kernel Library (MKL) on a system.

## Overview

`FindMKL.cmake` is a standard CMake "Find Module". It searches the filesystem for MKL header files and libraries, then exposes the results through standard CMake variables (`MKL_FOUND`, `MKL_INCLUDE_DIRS`, `MKL_LIBRARIES`).

## Step-by-Step Search Process

### 1. Determine Threading Model

Before searching for files, the module decides whether to look for the sequential or parallel version of MKL:

```cmake
option(MKL_PARALLEL "if mkl should be parallel" OFF)
```

| Option | Library Searched | Extra Dependencies |
|--------|------------------|--------------------|
| `MKL_PARALLEL=OFF` (default) | `mkl_sequential` | None |
| `MKL_PARALLEL=ON` | `mkl_intel_thread` | `iomp5` + `Threads` |

When parallel mode is enabled, the module also calls:
```cmake
find_package(Threads)
```
to locate the system threading library (usually `libpthread`).

### 2. Search for the Header File `mkl.h`

The module searches for the MKL header directory in two passes:

**Pass 1 — Hints-only search:**
It looks in directories pointed to by these CMake variables and environment variables (in order of priority):

| CMake Variable | Environment Variable |
|----------------|----------------------|
| `MKLROOT` | `$ENV{MKLROOT}` |
| `MKL_PATH` | `$ENV{MKL_PATH}` |
| `MKL_ROOT` | `$ENV{MKL_ROOT}` |

```cmake
find_path(MKL_INCLUDE_DIR mkl.h
    PATHS ${MKLROOT} ${MKL_PATH} ${MKL_ROOT}
          $ENV{MKLROOT} $ENV{MKL_PATH} $ENV{MKL_ROOT}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH)
```

**Pass 2 — System-wide fallback:**
If the hint search fails, it falls back to CMake's default system search paths:

```cmake
find_path(MKL_INCLUDE_DIR mkl.h PATH_SUFFIXES include)
```

If `mkl.h` is found, `MKL_INCLUDE_DIRS` is set to that directory.

### 3. Determine Architecture

Based on the host processor, the module sets the library subdirectory and naming suffix:

```cmake
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    set(__pathsuffix "lib/intel64")
    set(__libsfx _lp64)
else()
    set(__pathsuffix "lib/ia32")
    set(__libsfx "")
endif()
```

| Architecture | Library Subdirectory | Library Suffix |
|--------------|---------------------|----------------|
| x86_64 | `lib/intel64` | `_lp64` |
| ia32 (32-bit x86) | `lib/ia32` | (none) |

### 4. Search for the Core Libraries

Using the same root hints (`MKLROOT`, `MKL_PATH`, etc.), the module searches for three core libraries:

**a) Main Intel MKL interface library:**
```cmake
find_library(MKL_LIB_INTEL
    PATHS ${MKLROOT} ...
    PATH_SUFFIXES lib ${__pathsuffix}
    NAMES mkl_intel${__libsfx})
```
Searches for: `libmkl_intel_lp64.so` (x86_64) or `libmkl_intel.so` (ia32)

**b) Threading/sequential layer:**
```cmake
find_library(MKL_LIB_SEQUENTIAL  ... NAMES mkl_sequential)
# or
find_library(MKL_LIB_INTEL_THREAD ... NAMES mkl_intel_thread)
```

**c) Computational core:**
```cmake
find_library(MKL_LIB_CORE ... NAMES mkl_core)
```
Searches for: `libmkl_core.so`

**d) OpenMP runtime (parallel mode only):**
```cmake
find_library(MKL_LIB_IOMP5 ... NAMES iomp5)
```
Searches for: `libiomp5.so`

### 5. Assemble the Final Library List

If all required core libraries are found, they are combined into `MKL_LIBRARIES`:

```cmake
set(MKL_LIBRARIES
    ${MKL_LIB_INTEL}
    ${MKL_LIB_SEQUENTIAL}      # or ${MKL_LIB_INTEL_THREAD}
    ${MKL_LIB_CORE}
    ${MKL_LIB_IOMP5}           # only if MKL_PARALLEL=ON
    ${CMAKE_THREAD_LIBS_INIT}  # only if MKL_PARALLEL=ON
)
```

### 6. Report Results

The module uses CMake's standard helper to validate and report:

```cmake
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG
    MKL_LIBRARIES MKL_INCLUDE_DIRS)
```

This sets `MKL_FOUND` to `TRUE` only if both `MKL_LIBRARIES` and `MKL_INCLUDE_DIRS` are non-empty.

## Summary of Variables

| Output Variable | Description |
|-----------------|-------------|
| `MKL_FOUND` | TRUE if MKL was found successfully |
| `MKL_INCLUDE_DIRS` | Directory containing `mkl.h` |
| `MKL_LIBRARIES` | List of full paths to MKL libraries to link |
| `MKL_INCLUDE_DIR` | Internal: path to `mkl.h` (mark_as_advanced) |
| `MKL_LIB_INTEL` | Internal: path to `mkl_intel[_lp64]` |
| `MKL_LIB_SEQUENTIAL` | Internal: path to `mkl_sequential` |
| `MKL_LIB_INTEL_THREAD` | Internal: path to `mkl_intel_thread` |
| `MKL_LIB_CORE` | Internal: path to `mkl_core` |
| `MKL_LIB_IOMP5` | Internal: path to `iomp5` |

## How to Use in a CMake Project

```cmake
find_package(MKL REQUIRED)

if(MKL_FOUND)
    include_directories(${MKL_INCLUDE_DIRS})
    target_link_libraries(my_target ${MKL_LIBRARIES})
endif()
```

Or with modern imported targets (if ecbuild creates them):
```cmake
find_package(MKL REQUIRED)
target_link_libraries(my_target MKL::MKL)
```

## Typical MKL Installation Layout (Linux)

```
/opt/intel/mkl/
├── include/
│   └── mkl.h
└── lib/intel64/
    ├── libmkl_intel_lp64.so
    ├── libmkl_sequential.so
    ├── libmkl_intel_thread.so
    ├── libmkl_core.so
    └── libiomp5.so
```

The `MKLROOT` environment variable is usually set to `/opt/intel/mkl` or similar by Intel's `mklvars.sh` setup script.
