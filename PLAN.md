# Plan: FFTW3-Compatible Wrapper for FFT992

## 1. Goal
Create a drop-in Fortran replacement for the FFTW3 library that exposes the standard FFTW Fortran interface but routes all calls to the existing `fft992` / `set99` implementation.

## 2. Scope
- **Transforms**: real-to-complex (`r2c`) and complex-to-real (`c2r`) only.
- **Ranks**: 1-D only (batch / “many” is supported, multi-dimensional plans will abort).
- **Thread safety**: plan creation is **not** thread-safe (no locks).
- **Error handling**: illegal vector sizes abort via `CALL ABORT()`.
- **In-place / out-of-place**: both supported.

## 3. Files
| File | Action | Purpose |
|------|--------|---------|
| `parkind1.F90` | Unchanged | Kind parameters (`JPIM`, `JPRB`) |
| `set99.F90` | Unchanged | Factorisation and trig tables |
| `fft992.F90` | Unchanged | Core FFT engine |
| **`fftw3.f90`** | **New** | Wrapper module (FFTW-compatible API) |
| `test_fftw3.f90` | **New** | Regression / round-trip tests |

## 4. Plan Registry
FFTW3 uses opaque `INTEGER*8` plan handles.  
The wrapper maintains a module-level fixed-size array of derived types:

```fortran
TYPE :: FFTW_PLAN_T
  INTEGER(KIND=JPIM) :: N, LOT
  INTEGER(KIND=JPIM) :: ISIGN          ! FFTW_FORWARD (-1) or FFTW_BACKWARD (+1)
  LOGICAL            :: INPLACE
  INTEGER(KIND=JPIM) :: INC, JUMP      ! input strides
  INTEGER(KIND=JPIM) :: OST, ODIS      ! output strides (for copy-back)
  REAL(KIND=JPRB), ALLOCATABLE :: TRIGS(:)
  INTEGER(KIND=JPIM) :: IFAX(10)
END TYPE
```

A simple free-list gives the next available handle. Destroying a plan deallocates `TRIGS` and frees the slot.

## 5. API Surface
```fortran
SUBROUTINE dfftw_plan_dft_r2c_1d(plan, n, in, out, flags)
SUBROUTINE dfftw_plan_dft_c2r_1d(plan, n, in, out, flags)
SUBROUTINE dfftw_plan_many_dft_r2c(plan, rank, n, howmany, &
     in, inembed, istride, idist, &
     out, onembed, ostride, odist, flags)
SUBROUTINE dfftw_plan_many_dft_c2r(plan, rank, n, howmany, &
     in, inembed, istride, idist, &
     out, onembed, ostride, odist, flags)

SUBROUTINE dfftw_execute_dft_r2c(plan, in, out)
SUBROUTINE dfftw_execute_dft_c2r(plan, in, out)
SUBROUTINE dfftw_destroy_plan(plan)
```

Constants (`FFTW_ESTIMATE`, `FFTW_MEASURE`, `FFTW_FORWARD`, `FFTW_BACKWARD`, …) are defined with the same integer values as FFTW3.

## 6. In-Place Detection
The wrapper compares the memory addresses of `in` and `out` using the `LOC` intrinsic:
```fortran
INPLACE = (LOC(in(1)) == LOC(out(1)))
```
This matches FFTW’s semantics (same pointer = in-place) and does **not** require the caller to declare arrays with `TARGET`.

## 7. Data Layout & Strides
`fft992` stores complex data as interleaved reals (`A(0), B(0), A(1), B(1), …`).  
A Fortran `COMPLEX(KIND=KIND(0.0D0))` array is memory-identical to a `REAL(8)` array of twice the length, so no data conversion is needed.

### Stride convention mapping (FFTW → wrapper)
- `r2c`: `in` is REAL, `out` is COMPLEX.
  - Input strides `istride`, `idist` are in **real** elements.
  - Output strides `ostride`, `odist` are in **complex** elements → internally multiplied by 2 when treated as a real buffer.
- `c2r`: the reverse.

### Execution paths
**In-place execution**: call `fft992` directly on the user’s array with the provided strides.

**Out-of-place execution**:
1. Allocate a contiguous temporary real buffer of size `(N+2) * LOT`.
2. Copy the user’s input into the temp buffer, respecting input strides.
3. Call `fft992` on the temp buffer with `INC=1`, `JUMP=N+2`.
4. For **R2C only**, scale the entire temp buffer by `N` (see §8).
5. Copy the temp buffer back to the user’s output, respecting output strides.

## 8. Scaling Convention
`fft992` normalizes its forward transform (`ISIGN=-1`) by `1/N`, whereas FFTW is unnormalized.  
To preserve FFTW semantics:
- **`r2c`**: after `fft992(..., ISIGN=-1)`, multiply every element of the result by `N`.
- **`c2r`**: call `fft992(..., ISIGN=+1)` with no extra scaling.

With this convention, a forward transform followed by a backward transform returns the original data multiplied by `N`, exactly as FFTW does.

## 9. Unsupported Sizes
During planning, `SET99` is called. If `N` contains factors other than 2, 3, 4, 5, 6, 8 (or more than one factor 8), `SET99` writes an error and returns. The wrapper detects this and calls `ABORT()` immediately.

## 10. Explicit Interfaces for `fft992`
Since `fft992` is a standalone subroutine (not in a module), the wrapper module will contain an `INTERFACE` block providing an explicit interface to `fft992` and `set99`, so no existing files are modified.

## 11. Test Program
A small driver (`test_fftw3.f90`) will be provided that:
- creates 1-D and batch plans,
- runs round-trip (`r2c` then `c2r`) tests,
- verifies that the result equals `N * original`,
- checks a few known analytical transforms.
