PROGRAM DUMMY
  IMPLICIT NONE

  INCLUDE 'fftw3.f03'

  !---------------------------------------------------------------------
  ! DGEMM parameters
  !---------------------------------------------------------------------
  INTEGER, PARAMETER :: M = 2, N = 2, K = 3
  REAL(8) :: A(M,K), B(K,N), C(M,N)
  REAL(8) :: ALPHA, BETA
  INTEGER :: I, J
  EXTERNAL :: DGEMM

  !---------------------------------------------------------------------
  ! FFTW parameters
  !---------------------------------------------------------------------
  INTEGER, PARAMETER :: NFFT = 30
  INTEGER, PARAMETER :: NFFT2 = 48
  REAL(8), PARAMETER :: PI = 4.0D0 * ATAN(1.0D0)
  REAL(4), PARAMETER :: PIS = 4.0E0 * ATAN(1.0E0)

  REAL(8) :: R1(NFFT+2), ORIG1(NFFT)
  COMPLEX(8) :: C1(NFFT/2+1)
  REAL(8) :: R2(NFFT2), ORIG2(NFFT2)
  COMPLEX(8) :: C2(NFFT2/2+1)

  REAL(4) :: RS1(NFFT+2), ORIGS1(NFFT)
  COMPLEX(4) :: CS1(NFFT/2+1)

  INTEGER(8) :: PLAN_FWD, PLAN_BWD
  REAL(8) :: ERR, T
  REAL(4) :: ERRS, TS
  LOGICAL :: ALL_OK

  ALL_OK = .TRUE.

  !=====================================================================
  ! Test 1: DGEMM
  !=====================================================================
  PRINT '(A)', '--- Test 1: DGEMM ---'

  A(1,1) = 1.0D0; A(1,2) = 2.0D0; A(1,3) = 3.0D0
  A(2,1) = 4.0D0; A(2,2) = 5.0D0; A(2,3) = 6.0D0

  B(1,1) = 1.0D0; B(1,2) = 2.0D0
  B(2,1) = 3.0D0; B(2,2) = 4.0D0
  B(3,1) = 5.0D0; B(3,2) = 6.0D0

  ALPHA = 1.0D0
  BETA  = 0.0D0
  C = 0.0D0

  CALL DGEMM('N', 'N', M, N, K, ALPHA, A, M, B, K, BETA, C, M)

  IF (ABS(C(1,1)-22.0D0) < 1.0D-12 .AND. &
      ABS(C(1,2)-28.0D0) < 1.0D-12 .AND. &
      ABS(C(2,1)-49.0D0) < 1.0D-12 .AND. &
      ABS(C(2,2)-64.0D0) < 1.0D-12) THEN
    PRINT '(A)', 'PASSED'
  ELSE
    PRINT '(A)', 'FAILED'
    ALL_OK = .FALSE.
  END IF

  !=====================================================================
  ! Test 2: DP FFTW in-place round-trip N=30
  !=====================================================================
  PRINT '(A)', '--- Test 2: DP FFTW in-place round-trip N=30 ---'
  DO I = 1, NFFT
    T = 2.0D0 * PI * DBLE(I-1) / DBLE(NFFT)
    ORIG1(I) = SIN(3.0D0 * T) + 0.5D0 * COS(5.0D0 * T)
    R1(I) = ORIG1(I)
  END DO

  CALL dfftw_plan_dft_r2c_1d(PLAN_FWD, NFFT, R1, C1, FFTW_ESTIMATE)
  CALL dfftw_execute_dft_r2c(PLAN_FWD, R1, C1)
  CALL dfftw_destroy_plan(PLAN_FWD)

  CALL dfftw_plan_dft_c2r_1d(PLAN_BWD, NFFT, C1, R1, FFTW_ESTIMATE)
  CALL dfftw_execute_dft_c2r(PLAN_BWD, C1, R1)
  CALL dfftw_destroy_plan(PLAN_BWD)

  ERR = 0.0D0
  DO I = 1, NFFT
    ERR = MAX(ERR, ABS(R1(I) - ORIG1(I) * DBLE(NFFT)))
  END DO
  PRINT '(A,ES12.4)', 'Max error = ', ERR
  IF (ERR > 1.0D-12) THEN
    PRINT '(A)', 'FAILED'
    ALL_OK = .FALSE.
  ELSE
    PRINT '(A)', 'PASSED'
  END IF

  !=====================================================================
  ! Test 3: DP FFTW out-of-place round-trip N=48
  !=====================================================================
  PRINT '(A)', '--- Test 3: DP FFTW out-of-place round-trip N=48 ---'
  DO I = 1, NFFT2
    T = 2.0D0 * PI * DBLE(I-1) / DBLE(NFFT2)
    ORIG2(I) = COS(2.0D0 * T) - 0.3D0 * SIN(7.0D0 * T)
    R2(I) = ORIG2(I)
  END DO

  CALL dfftw_plan_dft_r2c_1d(PLAN_FWD, NFFT2, R2, C2, FFTW_ESTIMATE)
  CALL dfftw_execute_dft_r2c(PLAN_FWD, R2, C2)
  CALL dfftw_destroy_plan(PLAN_FWD)

  CALL dfftw_plan_dft_c2r_1d(PLAN_BWD, NFFT2, C2, R2, FFTW_ESTIMATE)
  CALL dfftw_execute_dft_c2r(PLAN_BWD, C2, R2)
  CALL dfftw_destroy_plan(PLAN_BWD)

  ERR = 0.0D0
  DO I = 1, NFFT2
    ERR = MAX(ERR, ABS(R2(I) - ORIG2(I) * DBLE(NFFT2)))
  END DO
  PRINT '(A,ES12.4)', 'Max error = ', ERR
  IF (ERR > 1.0D-12) THEN
    PRINT '(A)', 'FAILED'
    ALL_OK = .FALSE.
  ELSE
    PRINT '(A)', 'PASSED'
  END IF

  !=====================================================================
  ! Test 4: SP FFTW round-trip N=30
  !=====================================================================
  PRINT '(A)', '--- Test 4: SP FFTW round-trip N=30 ---'
  DO I = 1, NFFT
    TS = 2.0E0 * PIS * REAL(I-1,4) / REAL(NFFT,4)
    ORIGS1(I) = SIN(3.0E0 * TS) + 0.5E0 * COS(5.0E0 * TS)
    RS1(I) = ORIGS1(I)
  END DO

  CALL sfftw_plan_dft_r2c_1d(PLAN_FWD, NFFT, RS1, CS1, FFTW_ESTIMATE)
  CALL sfftw_execute_dft_r2c(PLAN_FWD, RS1, CS1)
  CALL sfftw_destroy_plan(PLAN_FWD)

  CALL sfftw_plan_dft_c2r_1d(PLAN_BWD, NFFT, CS1, RS1, FFTW_ESTIMATE)
  CALL sfftw_execute_dft_c2r(PLAN_BWD, CS1, RS1)
  CALL sfftw_destroy_plan(PLAN_BWD)

  ERRS = 0.0E0
  DO I = 1, NFFT
    ERRS = MAX(ERRS, ABS(RS1(I) - ORIGS1(I) * REAL(NFFT,4)))
  END DO
  PRINT '(A,ES12.4)', 'Max error = ', DBLE(ERRS)
  IF (ERRS > 1.0E-5) THEN
    PRINT '(A)', 'FAILED'
    ALL_OK = .FALSE.
  ELSE
    PRINT '(A)', 'PASSED'
  END IF

  !=====================================================================
  ! Summary
  !=====================================================================
  IF (ALL_OK) THEN
    PRINT '(A)', '=== ALL TESTS PASSED ==='
    STOP 0
  ELSE
    PRINT '(A)', '=== SOME TESTS FAILED ==='
    STOP 1
  END IF

END PROGRAM DUMMY
