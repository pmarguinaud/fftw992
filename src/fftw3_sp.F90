!=======================================================================
! FFTW3 single-precision Fortran wrapper (libfftw3f)
!=======================================================================

MODULE FFTW3_SP_M
  IMPLICIT NONE

  ! Complex kind
  INTEGER, PARAMETER :: JPCM = KIND((1.0_4, 0.0_4))

  ! FFTW constants
  INTEGER(KIND=8), PARAMETER :: FFTW_FORWARD    = -1_8
  INTEGER(KIND=8), PARAMETER :: FFTW_BACKWARD   = +1_8
  INTEGER(KIND=8), PARAMETER :: FFTW_ESTIMATE   =  0_8
  INTEGER(KIND=8), PARAMETER :: FFTW_MEASURE    =  0_8
  INTEGER(KIND=8), PARAMETER :: FFTW_PATIENT    =  0_8
  INTEGER(KIND=8), PARAMETER :: FFTW_EXHAUSTIVE =  0_8

  ! Plan registry
  INTEGER, PARAMETER :: MAXPLANS = 1024
  TYPE :: PLAN_T
    LOGICAL            :: ACTIVE = .FALSE.
    INTEGER(KIND=4) :: N      = 0
    INTEGER(KIND=4) :: LOT    = 0
    INTEGER(KIND=4) :: ISIGN  = 0
    LOGICAL            :: INPLACE= .FALSE.
    INTEGER(KIND=4) :: INC    = 0
    INTEGER(KIND=4) :: JUMP   = 0
    INTEGER(KIND=4) :: OINC   = 0
    INTEGER(KIND=4) :: OJUMP  = 0
    REAL(KIND=4), ALLOCATABLE :: TRIGS(:)
    INTEGER(KIND=4) :: IFAX(10) = 0
  END TYPE
  TYPE(PLAN_T), SAVE :: PLANS(MAXPLANS)

  ! C library interfaces
  INTERFACE
    FUNCTION c_malloc(n) BIND(C, NAME='malloc')
      USE ISO_C_BINDING, ONLY: C_PTR, C_SIZE_T
      IMPLICIT NONE
      INTEGER(C_SIZE_T), VALUE :: n
      TYPE(C_PTR) :: c_malloc
    END FUNCTION c_malloc

    SUBROUTINE c_free(ptr) BIND(C, NAME='free')
      USE ISO_C_BINDING, ONLY: C_PTR
      IMPLICIT NONE
      TYPE(C_PTR), VALUE :: ptr
    END SUBROUTINE c_free
  END INTERFACE

  ! Internal interfaces to the core FFT routines
  INTERFACE
    SUBROUTINE FFT992_SP(A,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
      REAL(KIND=4), INTENT(INOUT) :: A(*)
      REAL(KIND=4), INTENT(IN)    :: TRIGS(N)
      INTEGER(KIND=4), INTENT(IN) :: IFAX(10)
      INTEGER(KIND=4), INTENT(IN) :: INC, JUMP, N, LOT, ISIGN
    END SUBROUTINE FFT992_SP

    SUBROUTINE SET99_SP(TRIGS,IFAX,N)
      REAL(KIND=4), INTENT(OUT)   :: TRIGS(N)
      INTEGER(KIND=4), INTENT(OUT):: IFAX(*)
      INTEGER(KIND=4), INTENT(IN) :: N
    END SUBROUTINE SET99_SP
  END INTERFACE

CONTAINS

  !---------------------------------------------------------------------
  ! In-place detection
  !---------------------------------------------------------------------
  LOGICAL FUNCTION IP_RC_SP(RIN,COUT)
    REAL(KIND=4), INTENT(IN) :: RIN(*)
    COMPLEX(KIND=JPCM), INTENT(IN) :: COUT(*)
    IP_RC_SP = (LOC(RIN(1)) == LOC(COUT(1)))
  END FUNCTION IP_RC_SP

  LOGICAL FUNCTION IP_CR_SP(CIN,ROUT)
    COMPLEX(KIND=JPCM), INTENT(IN) :: CIN(*)
    REAL(KIND=4), INTENT(IN) :: ROUT(*)
    IP_CR_SP = (LOC(CIN(1)) == LOC(ROUT(1)))
  END FUNCTION IP_CR_SP

  !---------------------------------------------------------------------
  ! Register a new plan
  !---------------------------------------------------------------------
  SUBROUTINE REG_PLAN(N,LOT,ISIGN,INPLACE,INC,JUMP,OINC,OJUMP, &
                      TRIGS,IFAX,HANDLE)
    INTEGER(KIND=4), INTENT(IN) :: N,LOT,ISIGN
    LOGICAL,         INTENT(IN) :: INPLACE
    INTEGER(KIND=4), INTENT(IN) :: INC,JUMP,OINC,OJUMP
    REAL(KIND=4),    INTENT(IN) :: TRIGS(N)
    INTEGER(KIND=4), INTENT(IN) :: IFAX(10)
    INTEGER(KIND=8), INTENT(OUT):: HANDLE

    INTEGER(KIND=8) :: I

    DO I = 1, MAXPLANS
      IF (.NOT. PLANS(I)%ACTIVE) THEN
        PLANS(I)%ACTIVE  = .TRUE.
        PLANS(I)%N       = N
        PLANS(I)%LOT     = LOT
        PLANS(I)%ISIGN   = ISIGN
        PLANS(I)%INPLACE = INPLACE
        PLANS(I)%INC     = INC
        PLANS(I)%JUMP    = JUMP
        PLANS(I)%OINC    = OINC
        PLANS(I)%OJUMP   = OJUMP
        PLANS(I)%IFAX    = IFAX
        ALLOCATE(PLANS(I)%TRIGS(N))
        PLANS(I)%TRIGS   = TRIGS
        HANDLE = I
        RETURN
      END IF
    END DO

    WRITE(0,*) 'FFTW3: plan registry full (max=',MAXPLANS,')'
    CALL ABORT()
  END SUBROUTINE REG_PLAN

  !---------------------------------------------------------------------
  ! Validate N and build trig tables
  !---------------------------------------------------------------------
  SUBROUTINE PREPARE_TRIGS(N,TRIGS,IFAX)
    INTEGER(KIND=4), INTENT(IN)  :: N
    REAL(KIND=4),    INTENT(OUT) :: TRIGS(N)
    INTEGER(KIND=4), INTENT(OUT) :: IFAX(10)

    IF (N < 2) THEN
      WRITE(0,*) 'FFTW3: N must be >= 2, got ',N
      CALL ABORT()
    END IF

    CALL SET99_SP(TRIGS,IFAX,N)

    IF (IFAX(1) <= 0) THEN
      WRITE(0,*) 'FFTW3: N =',N,' contains illegal factors for fft992'
      CALL ABORT()
    END IF
  END SUBROUTINE PREPARE_TRIGS

  !---------------------------------------------------------------------
  ! Scale in-place R2C result by N
  !---------------------------------------------------------------------
  SUBROUTINE SCALE_R2C_SP(A,N,LOT,INC,JUMP)
    REAL(KIND=4), INTENT(INOUT) :: A(*)
    INTEGER(KIND=4), INTENT(IN) :: N,LOT,INC,JUMP
    INTEGER(KIND=4) :: J,K,II
    REAL(KIND=4) :: SCALE
    SCALE = REAL(N,KIND=4)
    DO J = 1, LOT
      II = 1 + (J-1)*JUMP
      DO K = 1, N+2
        A(II) = A(II)*SCALE
        II = II + INC
      END DO
    END DO
  END SUBROUTINE SCALE_R2C_SP

END MODULE FFTW3_SP_M


!=======================================================================
! Public API subroutines for libfftw3f
!=======================================================================

SUBROUTINE sfftw_plan_dft_r2c_1d(PLAN,N,IN,OUT,FLAGS)
  USE FFTW3_SP_M
  INTEGER(KIND=8), INTENT(OUT) :: PLAN
  INTEGER,         INTENT(IN)  :: N
  REAL(KIND=4),    INTENT(INOUT) :: IN(*)
  COMPLEX(KIND=JPCM), INTENT(INOUT) :: OUT(*)
  INTEGER(KIND=8), INTENT(IN)  :: FLAGS
  INTEGER(KIND=4) :: IFAX(10)
  REAL(KIND=4)    :: TRIGS(N)
  LOGICAL         :: INPLACE
  NN = INT(N,KIND=4)
  CALL PREPARE_TRIGS(NN,TRIGS,IFAX)
  INPLACE = IP_RC_SP(IN,OUT)
  CALL REG_PLAN(NN,1_4,-1_4,INPLACE,1_4,NN+2,2_4,NN+2, &
                TRIGS,IFAX,PLAN)
END SUBROUTINE sfftw_plan_dft_r2c_1d

SUBROUTINE sfftw_plan_dft_c2r_1d(PLAN,N,IN,OUT,FLAGS)
  USE FFTW3_SP_M
  INTEGER(KIND=8), INTENT(OUT) :: PLAN
  INTEGER,         INTENT(IN)  :: N
  COMPLEX(KIND=JPCM), INTENT(INOUT) :: IN(*)
  REAL(KIND=4),    INTENT(INOUT) :: OUT(*)
  INTEGER(KIND=8), INTENT(IN)  :: FLAGS
  INTEGER(KIND=4) :: IFAX(10)
  REAL(KIND=4)    :: TRIGS(N)
  LOGICAL         :: INPLACE
  NN = INT(N,KIND=4)
  CALL PREPARE_TRIGS(NN,TRIGS,IFAX)
  INPLACE = IP_CR_SP(IN,OUT)
  CALL REG_PLAN(NN,1_4,+1_4,INPLACE,1_4,NN+2,1_4,NN+2, &
                TRIGS,IFAX,PLAN)
END SUBROUTINE sfftw_plan_dft_c2r_1d

SUBROUTINE sfftw_plan_many_dft_r2c(PLAN,RANK,N,HOWMANY, &
     IN,INEMBED,ISTRIDE,IDIST, &
     OUT,ONEMBED,OSTRIDE,ODIST,FLAGS)
  USE FFTW3_SP_M
  INTEGER(KIND=8), INTENT(OUT) :: PLAN
  INTEGER,         INTENT(IN)  :: RANK,N(*),HOWMANY
  REAL(KIND=4),    INTENT(INOUT) :: IN(*)
  INTEGER,         INTENT(IN)  :: INEMBED(*),ISTRIDE,IDIST
  COMPLEX(KIND=JPCM), INTENT(INOUT) :: OUT(*)
  INTEGER,         INTENT(IN)  :: ONEMBED(*),OSTRIDE,ODIST
  INTEGER(KIND=8), INTENT(IN)  :: FLAGS
  INTEGER(KIND=4) :: IFAX(10)
  REAL(KIND=4)    :: TRIGS(N(1))
  LOGICAL         :: INPLACE
  IF (RANK /= 1) THEN; WRITE(0,*) 'FFTW3: only rank=1'; CALL ABORT(); END IF
  CALL PREPARE_TRIGS(INT(N(1),4),TRIGS,IFAX)
  INPLACE = IP_RC_SP(IN,OUT)
  CALL REG_PLAN(INT(N(1),4),INT(HOWMANY,4),-1_4,INPLACE, &
                INT(ISTRIDE,4),INT(IDIST,4), &
                INT(OSTRIDE,4)*2,INT(ODIST,4)*2, &
                TRIGS,IFAX,PLAN)
END SUBROUTINE sfftw_plan_many_dft_r2c

SUBROUTINE sfftw_plan_many_dft_c2r(PLAN,RANK,N,HOWMANY, &
     IN,INEMBED,ISTRIDE,IDIST, &
     OUT,ONEMBED,OSTRIDE,ODIST,FLAGS)
  USE FFTW3_SP_M
  INTEGER(KIND=8), INTENT(OUT) :: PLAN
  INTEGER,         INTENT(IN)  :: RANK,N(*),HOWMANY
  COMPLEX(KIND=JPCM), INTENT(INOUT) :: IN(*)
  INTEGER,         INTENT(IN)  :: INEMBED(*),ISTRIDE,IDIST
  REAL(KIND=4),    INTENT(INOUT) :: OUT(*)
  INTEGER,         INTENT(IN)  :: ONEMBED(*),OSTRIDE,ODIST
  INTEGER(KIND=8), INTENT(IN)  :: FLAGS
  INTEGER(KIND=4) :: IFAX(10)
  REAL(KIND=4)    :: TRIGS(N(1))
  LOGICAL         :: INPLACE
  IF (RANK /= 1) THEN; WRITE(0,*) 'FFTW3: only rank=1'; CALL ABORT(); END IF
  CALL PREPARE_TRIGS(INT(N(1),4),TRIGS,IFAX)
  INPLACE = IP_CR_SP(IN,OUT)
  CALL REG_PLAN(INT(N(1),4),INT(HOWMANY,4),+1_4,INPLACE, &
                INT(ISTRIDE,4)*2,INT(IDIST,4)*2, &
                INT(OSTRIDE,4),INT(ODIST,4), &
                TRIGS,IFAX,PLAN)
END SUBROUTINE sfftw_plan_many_dft_c2r

SUBROUTINE sfftw_execute_dft_r2c(PLAN,IN,OUT)
  USE FFTW3_SP_M
  INTEGER(KIND=8), INTENT(IN) :: PLAN
  REAL(KIND=4),    INTENT(INOUT) :: IN(*)
  COMPLEX(KIND=JPCM), INTENT(INOUT) :: OUT(*)
  INTEGER(KIND=8)    :: ID
  INTEGER(KIND=4) :: N,LOT
  REAL(KIND=4), ALLOCATABLE :: WORK(:)
  INTEGER(KIND=4) :: J,K,II,IO
  ID = PLAN
  IF (ID < 1 .OR. ID > MAXPLANS .OR. .NOT. PLANS(ID)%ACTIVE) THEN
    WRITE(0,*) 'FFTW3: invalid SP plan in execute_r2c'; CALL ABORT()
  END IF
  N   = PLANS(ID)%N
  LOT = PLANS(ID)%LOT
  IF (PLANS(ID)%INPLACE) THEN
    CALL FFT992_SP(IN,PLANS(ID)%TRIGS,PLANS(ID)%IFAX, &
                   PLANS(ID)%INC,PLANS(ID)%JUMP,N,LOT,-1_4)
    CALL SCALE_R2C_SP(IN,N,LOT,PLANS(ID)%INC,PLANS(ID)%JUMP)
  ELSE
    ALLOCATE(WORK((N+2)*LOT))
    DO J = 1, LOT
      II = 1 + (J-1)*PLANS(ID)%JUMP
      IO = 1 + (J-1)*(N+2)
      DO K = 1, N
        WORK(IO) = IN(II); II = II + PLANS(ID)%INC; IO = IO + 1
      END DO
      WORK(IO) = 0.0_4; WORK(IO+1) = 0.0_4
    END DO
    CALL FFT992_SP(WORK,PLANS(ID)%TRIGS,PLANS(ID)%IFAX, &
                   1_4,N+2,N,LOT,-1_4)
    DO J = 1, (N+2)*LOT
      WORK(J) = WORK(J)*REAL(N,KIND=4)
    END DO
    DO J = 1, LOT
      II = 1 + (J-1)*(N+2)
      IO = 1 + (J-1)*(PLANS(ID)%OJUMP/2)
      DO K = 1, N/2 + 1
        OUT(IO) = CMPLX(WORK(II),WORK(II+1),KIND=JPCM)
        II = II + 2; IO = IO + PLANS(ID)%OINC/2
      END DO
    END DO
    DEALLOCATE(WORK)
  END IF
END SUBROUTINE sfftw_execute_dft_r2c

SUBROUTINE sfftw_execute_dft_c2r(PLAN,IN,OUT)
  USE FFTW3_SP_M
  INTEGER(KIND=8), INTENT(IN) :: PLAN
  COMPLEX(KIND=JPCM), INTENT(INOUT) :: IN(*)
  REAL(KIND=4),    INTENT(INOUT) :: OUT(*)
  INTEGER(KIND=8)    :: ID
  INTEGER(KIND=4) :: N,LOT
  REAL(KIND=4), ALLOCATABLE :: WORK(:)
  INTEGER(KIND=4) :: J,K,II,IO
  ID = PLAN
  IF (ID < 1 .OR. ID > MAXPLANS .OR. .NOT. PLANS(ID)%ACTIVE) THEN
    WRITE(0,*) 'FFTW3: invalid SP plan in execute_c2r'; CALL ABORT()
  END IF
  N   = PLANS(ID)%N
  LOT = PLANS(ID)%LOT
  IF (PLANS(ID)%INPLACE) THEN
    CALL FFT992_SP(OUT,PLANS(ID)%TRIGS,PLANS(ID)%IFAX, &
                   1_4,N+2,N,LOT,+1_4)
  ELSE
    ALLOCATE(WORK((N+2)*LOT))
    DO J = 1, LOT
      II = 1 + (J-1)*(PLANS(ID)%JUMP/2)
      IO = 1 + (J-1)*(N+2)
      DO K = 1, N/2 + 1
        WORK(IO)   = REAL(IN(II),KIND=4)
        WORK(IO+1) = AIMAG(IN(II))
        II = II + MAX(PLANS(ID)%INC/2,1_4)
        IO = IO + 2
      END DO
    END DO
    CALL FFT992_SP(WORK,PLANS(ID)%TRIGS,PLANS(ID)%IFAX, &
                   1_4,N+2,N,LOT,+1_4)
    DO J = 1, LOT
      II = 1 + (J-1)*(N+2)
      IO = 1 + (J-1)*PLANS(ID)%OJUMP
      DO K = 1, N
        OUT(IO) = WORK(II); II = II + 1; IO = IO + PLANS(ID)%OINC
      END DO
    END DO
    DEALLOCATE(WORK)
  END IF
END SUBROUTINE sfftw_execute_dft_c2r

SUBROUTINE sfftw_destroy_plan(PLAN)
  USE FFTW3_SP_M
  INTEGER(KIND=8), INTENT(IN) :: PLAN
  INTEGER(KIND=8) :: ID
  ID = PLAN
  IF (ID >= 1 .AND. ID <= MAXPLANS) THEN
    IF (ALLOCATED(PLANS(ID)%TRIGS)) DEALLOCATE(PLANS(ID)%TRIGS)
    PLANS(ID)%ACTIVE = .FALSE.
  END IF
END SUBROUTINE sfftw_destroy_plan

FUNCTION fftwf_alloc_complex(N) RESULT(PTR)
  USE ISO_C_BINDING, ONLY: C_PTR, C_SIZE_T
  USE FFTW3_SP_M, ONLY: c_malloc
  INTEGER(C_SIZE_T), INTENT(IN) :: N
  TYPE(C_PTR) :: PTR
  PTR = c_malloc((N + 1) * 8_C_SIZE_T)
END FUNCTION fftwf_alloc_complex

SUBROUTINE fftwf_free(P)
  USE ISO_C_BINDING, ONLY: C_PTR
  USE FFTW3_SP_M, ONLY: c_free
  TYPE(C_PTR), INTENT(IN) :: P
  CALL c_free(P)
END SUBROUTINE fftwf_free
