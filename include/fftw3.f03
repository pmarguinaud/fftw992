!=======================================================================
! Public include file: FFTW3 user-facing interfaces and constants.
! Include this file in your program (or USE the FFTW3 module) to get
! explicit interfaces for the full DP + SP API.
!=======================================================================

  ! Complex kinds (self-contained, do not depend on PARKIND1)
  INTEGER, PARAMETER :: JPCB = KIND((1.0D0, 0.0D0))
  INTEGER, PARAMETER :: JPCM = KIND((1.0,  0.0))

  ! FFTW constants (same values as the C library)
  INTEGER(KIND=8), PARAMETER :: FFTW_FORWARD    = -1_8
  INTEGER(KIND=8), PARAMETER :: FFTW_BACKWARD   = +1_8
  INTEGER(KIND=8), PARAMETER :: FFTW_ESTIMATE   =  0_8
  INTEGER(KIND=8), PARAMETER :: FFTW_MEASURE    =  0_8
  INTEGER(KIND=8), PARAMETER :: FFTW_PATIENT    =  0_8
  INTEGER(KIND=8), PARAMETER :: FFTW_EXHAUSTIVE =  0_8
  INTEGER(KIND=8), PARAMETER :: FFTW_NO_SIMD    =  0_8

  !---------------------------------------------------------------------
  ! Memory allocation (ISO_C_BINDING types, matching FFTW3 F2003)
  !---------------------------------------------------------------------
  INTERFACE
    FUNCTION fftw_alloc_complex(n) RESULT(ptr)
      USE ISO_C_BINDING, ONLY: C_PTR, C_SIZE_T
      INTEGER(C_SIZE_T), INTENT(IN) :: n
      TYPE(C_PTR) :: ptr
    END FUNCTION fftw_alloc_complex

    SUBROUTINE fftw_free(p)
      USE ISO_C_BINDING, ONLY: C_PTR
      TYPE(C_PTR), INTENT(IN) :: p
    END SUBROUTINE fftw_free

    FUNCTION fftwf_alloc_complex(n) RESULT(ptr)
      USE ISO_C_BINDING, ONLY: C_PTR, C_SIZE_T
      INTEGER(C_SIZE_T), INTENT(IN) :: n
      TYPE(C_PTR) :: ptr
    END FUNCTION fftwf_alloc_complex

    SUBROUTINE fftwf_free(p)
      USE ISO_C_BINDING, ONLY: C_PTR
      TYPE(C_PTR), INTENT(IN) :: p
    END SUBROUTINE fftwf_free
  END INTERFACE
