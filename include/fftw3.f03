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

  !---------------------------------------------------------------------
  ! Double-precision API
  !---------------------------------------------------------------------
  INTERFACE
    SUBROUTINE dfftw_plan_dft_r2c_1d(plan, n, in, out, flags)
      INTEGER(KIND=8), INTENT(OUT) :: plan
      INTEGER,         INTENT(IN)  :: n
      REAL(KIND=8),    INTENT(INOUT) :: in(*)
      COMPLEX(KIND=8), INTENT(INOUT) :: out(*)
      INTEGER(KIND=8), INTENT(IN)  :: flags
    END SUBROUTINE dfftw_plan_dft_r2c_1d

    SUBROUTINE dfftw_plan_dft_c2r_1d(plan, n, in, out, flags)
      INTEGER(KIND=8), INTENT(OUT) :: plan
      INTEGER,         INTENT(IN)  :: n
      COMPLEX(KIND=8), INTENT(INOUT) :: in(*)
      REAL(KIND=8),    INTENT(INOUT) :: out(*)
      INTEGER(KIND=8), INTENT(IN)  :: flags
    END SUBROUTINE dfftw_plan_dft_c2r_1d

    SUBROUTINE dfftw_plan_many_dft_r2c(plan, rank, n, howmany, &
         in, inembed, istride, idist, &
         out, onembed, ostride, odist, flags)
      INTEGER(KIND=8), INTENT(OUT) :: plan
      INTEGER,         INTENT(IN)  :: rank
      INTEGER,         INTENT(IN)  :: n(*)
      INTEGER,         INTENT(IN)  :: howmany
      REAL(KIND=8),    INTENT(INOUT) :: in(*)
      INTEGER,         INTENT(IN)  :: inembed(*)
      INTEGER,         INTENT(IN)  :: istride
      INTEGER,         INTENT(IN)  :: idist
      COMPLEX(KIND=8), INTENT(INOUT) :: out(*)
      INTEGER,         INTENT(IN)  :: onembed(*)
      INTEGER,         INTENT(IN)  :: ostride
      INTEGER,         INTENT(IN)  :: odist
      INTEGER(KIND=8), INTENT(IN)  :: flags
    END SUBROUTINE dfftw_plan_many_dft_r2c

    SUBROUTINE dfftw_plan_many_dft_c2r(plan, rank, n, howmany, &
         in, inembed, istride, idist, &
         out, onembed, ostride, odist, flags)
      INTEGER(KIND=8), INTENT(OUT) :: plan
      INTEGER,         INTENT(IN)  :: rank
      INTEGER,         INTENT(IN)  :: n(*)
      INTEGER,         INTENT(IN)  :: howmany
      COMPLEX(KIND=8), INTENT(INOUT) :: in(*)
      INTEGER,         INTENT(IN)  :: inembed(*)
      INTEGER,         INTENT(IN)  :: istride
      INTEGER,         INTENT(IN)  :: idist
      REAL(KIND=8),    INTENT(INOUT) :: out(*)
      INTEGER,         INTENT(IN)  :: onembed(*)
      INTEGER,         INTENT(IN)  :: ostride
      INTEGER,         INTENT(IN)  :: odist
      INTEGER(KIND=8), INTENT(IN)  :: flags
    END SUBROUTINE dfftw_plan_many_dft_c2r

    SUBROUTINE dfftw_execute_dft_r2c(plan, in, out)
      INTEGER(KIND=8), INTENT(IN) :: plan
      REAL(KIND=8),    INTENT(INOUT) :: in(*)
      COMPLEX(KIND=8), INTENT(INOUT) :: out(*)
    END SUBROUTINE dfftw_execute_dft_r2c

    SUBROUTINE dfftw_execute_dft_c2r(plan, in, out)
      INTEGER(KIND=8), INTENT(IN) :: plan
      COMPLEX(KIND=8), INTENT(INOUT) :: in(*)
      REAL(KIND=8),    INTENT(INOUT) :: out(*)
    END SUBROUTINE dfftw_execute_dft_c2r

    SUBROUTINE dfftw_destroy_plan(plan)
      INTEGER(KIND=8), INTENT(IN) :: plan
    END SUBROUTINE dfftw_destroy_plan
  END INTERFACE

  !---------------------------------------------------------------------
  ! Single-precision API
  !---------------------------------------------------------------------
  INTERFACE
    SUBROUTINE sfftw_plan_dft_r2c_1d(plan, n, in, out, flags)
      INTEGER(KIND=8), INTENT(OUT) :: plan
      INTEGER,         INTENT(IN)  :: n
      REAL(KIND=4),    INTENT(INOUT) :: in(*)
      COMPLEX(KIND=4), INTENT(INOUT) :: out(*)
      INTEGER(KIND=8), INTENT(IN)  :: flags
    END SUBROUTINE sfftw_plan_dft_r2c_1d

    SUBROUTINE sfftw_plan_dft_c2r_1d(plan, n, in, out, flags)
      INTEGER(KIND=8), INTENT(OUT) :: plan
      INTEGER,         INTENT(IN)  :: n
      COMPLEX(KIND=4), INTENT(INOUT) :: in(*)
      REAL(KIND=4),    INTENT(INOUT) :: out(*)
      INTEGER(KIND=8), INTENT(IN)  :: flags
    END SUBROUTINE sfftw_plan_dft_c2r_1d

    SUBROUTINE sfftw_plan_many_dft_r2c(plan, rank, n, howmany, &
         in, inembed, istride, idist, &
         out, onembed, ostride, odist, flags)
      INTEGER(KIND=8), INTENT(OUT) :: plan
      INTEGER,         INTENT(IN)  :: rank
      INTEGER,         INTENT(IN)  :: n(*)
      INTEGER,         INTENT(IN)  :: howmany
      REAL(KIND=4),    INTENT(INOUT) :: in(*)
      INTEGER,         INTENT(IN)  :: inembed(*)
      INTEGER,         INTENT(IN)  :: istride
      INTEGER,         INTENT(IN)  :: idist
      COMPLEX(KIND=4), INTENT(INOUT) :: out(*)
      INTEGER,         INTENT(IN)  :: onembed(*)
      INTEGER,         INTENT(IN)  :: ostride
      INTEGER,         INTENT(IN)  :: odist
      INTEGER(KIND=8), INTENT(IN)  :: flags
    END SUBROUTINE sfftw_plan_many_dft_r2c

    SUBROUTINE sfftw_plan_many_dft_c2r(plan, rank, n, howmany, &
         in, inembed, istride, idist, &
         out, onembed, ostride, odist, flags)
      INTEGER(KIND=8), INTENT(OUT) :: plan
      INTEGER,         INTENT(IN)  :: rank
      INTEGER,         INTENT(IN)  :: n(*)
      INTEGER,         INTENT(IN)  :: howmany
      COMPLEX(KIND=4), INTENT(INOUT) :: in(*)
      INTEGER,         INTENT(IN)  :: inembed(*)
      INTEGER,         INTENT(IN)  :: istride
      INTEGER,         INTENT(IN)  :: idist
      REAL(KIND=4),    INTENT(INOUT) :: out(*)
      INTEGER,         INTENT(IN)  :: onembed(*)
      INTEGER,         INTENT(IN)  :: ostride
      INTEGER,         INTENT(IN)  :: odist
      INTEGER(KIND=8), INTENT(IN)  :: flags
    END SUBROUTINE sfftw_plan_many_dft_c2r

    SUBROUTINE sfftw_execute_dft_r2c(plan, in, out)
      INTEGER(KIND=8), INTENT(IN) :: plan
      REAL(KIND=4),    INTENT(INOUT) :: in(*)
      COMPLEX(KIND=4), INTENT(INOUT) :: out(*)
    END SUBROUTINE sfftw_execute_dft_r2c

    SUBROUTINE sfftw_execute_dft_c2r(plan, in, out)
      INTEGER(KIND=8), INTENT(IN) :: plan
      COMPLEX(KIND=4), INTENT(INOUT) :: in(*)
      REAL(KIND=4),    INTENT(INOUT) :: out(*)
    END SUBROUTINE sfftw_execute_dft_c2r

    SUBROUTINE sfftw_destroy_plan(plan)
      INTEGER(KIND=8), INTENT(IN) :: plan
    END SUBROUTINE sfftw_destroy_plan
  END INTERFACE

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
