# How ecbuild Finds FFTW

This document summarises the search logic in `ecbuild/cmake/FindFFTW.cmake`.

## Entry Point

```cmake
find_package(FFTW [REQUIRED] [QUIET]
             [COMPONENTS [single] [double] [long_double] [quad]])
```

Default component (if none specified): `double`.

---

## Search Priority (tried in order)

### 1. User override
If the caller has already set:
- `FFTW_LIBRARIES`
- `FFTW_INCLUDE_DIRS`

…then **no searching is performed** and those values are used directly.

---

### 2. MKL implementation (preferred by default)

**Default**: MKL is tried first unless explicitly disabled.

| Condition | Behaviour |
|-----------|-----------|
| `FFTW_ENABLE_MKL=ON` (explicit) | Only MKL is considered |
| `FFTW_ENABLE_MKL=OFF` (explicit) | MKL is skipped |
| Neither defined | MKL is preferred unless `ENABLE_MKL=OFF` |

If MKL is enabled, it calls:
```cmake
find_package(MKL)
```

If MKL is found:
- `FFTW_INCLUDE_DIRS = ${MKL_INCLUDE_DIRS}/fftw`
- `FFTW_LIBRARIES     = ${MKL_LIBRARIES}`
- All requested components are marked found

The `MKLROOT` environment variable (or CMake vars `MKLROOT`, `MKL_PATH`, `MKL_ROOT`) is used by `FindMKL.cmake` to locate `mkl.h` and the MKL libraries.

---

### 3. ARMPL implementation

Same logic as MKL but for ARM Performance Libraries.

| Variable | Purpose |
|----------|---------|
| `FFTW_ENABLE_ARMPL` | Explicit ON/OFF |
| `ENABLE_ARMPL` | Fallback if `FFTW_ENABLE_ARMPL` undefined |

If ARMPL is found, `FFTW_INCLUDE_DIRS` and `FFTW_LIBRARIES` are set from `ARMPL_INCLUDE_DIRS` / `ARMPL_LIBRARIES`.

---

### 4. NVPL implementation

Same logic as ARMPL but for NVIDIA Performance Libraries.

| Variable | Purpose |
|----------|---------|
| `FFTW_ENABLE_NVPL` | Explicit ON/OFF |
| `ENABLE_NVPL` | Fallback if `FFTW_ENABLE_NVPL` undefined |

---

### 5. Standard FFTW library (final fallback)

Only reached if none of the vendor implementations above succeeded.

#### Search paths (in order)
1. `FFTW_ROOT` CMake variable
2. Environment variables: `FFTW_ROOT`, `FFTW_DIR`, `FFTWDIR`, `FFTW_PATH`
3. `FFTW_DIR`, `FFTW_PATH` CMake variables
4. `pkg-config` (`fftw3`) — if no `FFTW_ROOT` is set

Note: if the resolved directory name matches `lib`, it is treated as the **lib** directory and the parent is used as `FFTW_ROOT`.

#### What is searched
- **Include header**: `fftw3.h` under `${FFTW_ROOT}/include`
- **Libraries** (per requested component):
  - `fftw3`  → double precision
  - `fftw3f` → single precision
  - `fftw3l` → long double
  - `fftw3q` → quad precision

Library search paths: `${FFTW_ROOT}/lib`, `${FFTW_ROOT}/lib64`, or `pkg-config` library dirs.

---

## Output Variables

| Variable | Description |
|----------|-------------|
| `FFTW_FOUND` | TRUE if FFTW was found |
| `FFTW_INCLUDE_DIRS` | Include directory path |
| `FFTW_LIBRARIES` | Full paths to libraries to link |
| `FFTW_double_LIBRARIES` / `FFTW_single_LIBRARIES` / … | Per-component library paths |
| `FFTW::fftw3` / `FFTW::fftw3f` / … | Imported interface targets (modern CMake) |

---

## Key Notes

- `FindFFTW.cmake` **delegates to `FindMKL.cmake`** for the MKL case. `FindMKL.cmake` searches for `mkl.h` and libraries named `mkl_intel_lp64`, `mkl_sequential`, `mkl_core` (plus `iomp5` if `MKL_PARALLEL=ON`).
- The `/fftw` suffix on `FFTW_INCLUDE_DIRS` when using MKL means Intel MKL places `fftw3.f03` in `<mklroot>/include/fftw/`, not directly in `include/`.
- `FindFFTW.cmake` does **not** search for `fftw3.f03` directly — it searches for `fftw3.h` in the standard FFTW path, and for MKL it appends `/fftw` to the MKL include directory.
