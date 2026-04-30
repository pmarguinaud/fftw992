PROGRAM TEST_FFTW3
  USE ISO_C_BINDING, ONLY: C_PTR, C_ASSOCIATED, C_SIZE_T
  IMPLICIT NONE

  INCLUDE 'fftw3.f03'

  INTEGER, PARAMETER :: CK  = JPCB   ! DP complex kind
  INTEGER, PARAMETER :: CKS = JPCM   ! SP complex kind

  !---------------------------------------------------------------------
  ! Test parameters
  !---------------------------------------------------------------------
  INTEGER, PARAMETER :: N1 = 30
  INTEGER, PARAMETER :: N2 = 48
  INTEGER, PARAMETER :: LOT = 4
  REAL(KIND=8), PARAMETER :: PI = 4.0_8*ATAN(1.0_8)
  REAL(KIND=4), PARAMETER :: PIS = 4.0_4*ATAN(1.0_4)

  !---------------------------------------------------------------------
  ! DP local arrays
  !---------------------------------------------------------------------
  REAL(KIND=8)    :: A1(N1+2), ORIG1(N1)
  COMPLEX(KIND=CK)   :: C1(N1/2+1)
  REAL(KIND=8)    :: R1(N1+2)
  EQUIVALENCE (C1, R1)

  REAL(KIND=8)    :: A2(N2), B2(N2), ORIG2(N2)
  COMPLEX(KIND=CK)   :: C2(N2/2+1)

  REAL(KIND=8)    :: ABATCH(N1,LOT)
  COMPLEX(KIND=CK)   :: CBATCH(N1/2+1,LOT)
  REAL(KIND=8)    :: ORIG_BATCH(N1,LOT)

  REAL(KIND=8)    :: ABIP(N1+2,LOT)
  COMPLEX(KIND=CK)   :: CBIP(N1/2+1,LOT)
  REAL(KIND=8)    :: RBIP(N1+2,LOT)
  EQUIVALENCE (ABIP, CBIP, RBIP)

  !---------------------------------------------------------------------
  ! SP local arrays
  !---------------------------------------------------------------------
  REAL(KIND=4)    :: AS1(N1+2), ORIGS1(N1)
  COMPLEX(KIND=CKS)  :: CS1(N1/2+1)
  REAL(KIND=4)    :: RS1(N1+2)
  EQUIVALENCE (CS1, RS1)

  REAL(KIND=4)    :: AS2(N2), BS2(N2), ORIGS2(N2)
  COMPLEX(KIND=CKS)  :: CS2(N2/2+1)

  REAL(KIND=4)    :: ASBATCH(N1,LOT)
  COMPLEX(KIND=CKS)  :: CSBATCH(N1/2+1,LOT)
  REAL(KIND=4)    :: ORIGS_BATCH(N1,LOT)

  REAL(KIND=4)    :: ASBIP(N1+2,LOT)
  COMPLEX(KIND=CKS)  :: CSBIP(N1/2+1,LOT)
  REAL(KIND=4)    :: RSBIP(N1+2,LOT)
  EQUIVALENCE (ASBIP, CSBIP, RSBIP)

  !---------------------------------------------------------------------
  ! Memory handles
  !---------------------------------------------------------------------
  TYPE(C_PTR) :: HDP, HSP

  INTEGER(KIND=8) :: PLAN_FWD, PLAN_BWD
  INTEGER :: I,J
  REAL(KIND=8) :: ERR, T
  REAL(KIND=4) :: ERRS, TS
  LOGICAL :: ALL_OK

  ALL_OK = .TRUE.

  !=====================================================================
  ! Test DP 1: in-place R2C + C2R round-trip  (N=30)
  !=====================================================================
  PRINT '(A)', '--- DP Test 1: in-place round-trip N=30 ---'
  DO I = 1, N1
    T = 2.0_8*PI*REAL(I-1,KIND=8)/REAL(N1,KIND=8)
    ORIG1(I) = SIN(3.0_8*T) + 0.5_8*COS(5.0_8*T)
    R1(I) = ORIG1(I)
  END DO

  CALL dfftw_plan_dft_r2c_1d(PLAN_FWD,N1,R1,C1,FFTW_ESTIMATE)
  CALL dfftw_execute_dft_r2c(PLAN_FWD,R1,C1)
  CALL dfftw_destroy_plan(PLAN_FWD)

  CALL dfftw_plan_dft_c2r_1d(PLAN_BWD,N1,C1,R1,FFTW_ESTIMATE)
  CALL dfftw_execute_dft_c2r(PLAN_BWD,C1,R1)
  CALL dfftw_destroy_plan(PLAN_BWD)

  ERR = 0.0_8
  DO I = 1, N1
    ERR = MAX(ERR,ABS(R1(I)-ORIG1(I)*REAL(N1,KIND=8)))
  END DO
  PRINT '(A,ES12.4)', 'Max error = ', ERR
  IF (ERR > 1.0E-12_8) THEN; PRINT '(A)', 'FAILED'; ALL_OK = .FALSE.
  ELSE; PRINT '(A)', 'PASSED'; END IF

  !=====================================================================
  ! Test DP 2: out-of-place R2C + C2R round-trip (N=48)
  !=====================================================================
  PRINT '(A)', '--- DP Test 2: out-of-place round-trip N=48 ---'
  DO I = 1, N2
    T = 2.0_8*PI*REAL(I-1,KIND=8)/REAL(N2,KIND=8)
    ORIG2(I) = COS(2.0_8*T) - 0.3_8*SIN(7.0_8*T)
    A2(I) = ORIG2(I)
  END DO

  CALL dfftw_plan_dft_r2c_1d(PLAN_FWD,N2,A2,C2,FFTW_ESTIMATE)
  CALL dfftw_execute_dft_r2c(PLAN_FWD,A2,C2)
  CALL dfftw_destroy_plan(PLAN_FWD)

  CALL dfftw_plan_dft_c2r_1d(PLAN_BWD,N2,C2,B2,FFTW_ESTIMATE)
  CALL dfftw_execute_dft_c2r(PLAN_BWD,C2,B2)
  CALL dfftw_destroy_plan(PLAN_BWD)

  ERR = 0.0_8
  DO I = 1, N2
    ERR = MAX(ERR,ABS(B2(I)-ORIG2(I)*REAL(N2,KIND=8)))
  END DO
  PRINT '(A,ES12.4)', 'Max error = ', ERR
  IF (ERR > 1.0E-12_8) THEN; PRINT '(A)', 'FAILED'; ALL_OK = .FALSE.
  ELSE; PRINT '(A)', 'PASSED'; END IF

  !=====================================================================
  ! Test DP 3: batch out-of-place (N=30, LOT=4)
  !=====================================================================
  PRINT '(A)', '--- DP Test 3: batch out-of-place N=30, LOT=4 ---'
  DO J = 1, LOT
    DO I = 1, N1
      T = 2.0_8*PI*REAL(I-1,KIND=8)/REAL(N1,KIND=8)
      ORIG_BATCH(I,J) = SIN(REAL(J,KIND=8)*T)
      ABATCH(I,J) = ORIG_BATCH(I,J)
    END DO
  END DO

  CALL dfftw_plan_many_dft_r2c(PLAN_FWD,1,[N1],LOT, &
       ABATCH,[N1],1,N1, &
       CBATCH,[N1/2+1],1,N1/2+1, &
       FFTW_ESTIMATE)
  CALL dfftw_execute_dft_r2c(PLAN_FWD,ABATCH,CBATCH)
  CALL dfftw_destroy_plan(PLAN_FWD)

  CALL dfftw_plan_many_dft_c2r(PLAN_BWD,1,[N1],LOT, &
       CBATCH,[N1/2+1],1,N1/2+1, &
       ABATCH,[N1],1,N1, &
       FFTW_ESTIMATE)
  CALL dfftw_execute_dft_c2r(PLAN_BWD,CBATCH,ABATCH)
  CALL dfftw_destroy_plan(PLAN_BWD)

  ERR = 0.0_8
  DO J = 1, LOT
    DO I = 1, N1
      ERR = MAX(ERR,ABS(ABATCH(I,J)-ORIG_BATCH(I,J)*REAL(N1,KIND=8)))
    END DO
  END DO
  PRINT '(A,ES12.4)', 'Max error = ', ERR
  IF (ERR > 1.0E-12_8) THEN; PRINT '(A)', 'FAILED'; ALL_OK = .FALSE.
  ELSE; PRINT '(A)', 'PASSED'; END IF

  !=====================================================================
  ! Test DP 4: batch in-place (N=30, LOT=4)
  !=====================================================================
  PRINT '(A)', '--- DP Test 4: batch in-place N=30, LOT=4 ---'
  DO J = 1, LOT
    DO I = 1, N1
      T = 2.0_8*PI*REAL(I-1,KIND=8)/REAL(N1,KIND=8)
      ORIG_BATCH(I,J) = COS(REAL(J,KIND=8)*T)
      ABIP(I,J) = ORIG_BATCH(I,J)
    END DO
  END DO

  CALL dfftw_plan_many_dft_r2c(PLAN_FWD,1,[N1],LOT, &
       ABIP,[N1+2],1,N1+2, &
       CBIP,[N1/2+1],1,N1/2+1, &
       FFTW_ESTIMATE)
  CALL dfftw_execute_dft_r2c(PLAN_FWD,ABIP,CBIP)
  CALL dfftw_destroy_plan(PLAN_FWD)

  CALL dfftw_plan_many_dft_c2r(PLAN_BWD,1,[N1],LOT, &
       CBIP,[N1/2+1],1,N1/2+1, &
       RBIP,[N1+2],1,N1+2, &
       FFTW_ESTIMATE)
  CALL dfftw_execute_dft_c2r(PLAN_BWD,CBIP,RBIP)
  CALL dfftw_destroy_plan(PLAN_BWD)

  ERR = 0.0_8
  DO J = 1, LOT
    DO I = 1, N1
      ERR = MAX(ERR,ABS(RBIP(I,J)-ORIG_BATCH(I,J)*REAL(N1,KIND=8)))
    END DO
  END DO
  PRINT '(A,ES12.4)', 'Max error = ', ERR
  IF (ERR > 1.0E-12_8) THEN; PRINT '(A)', 'FAILED'; ALL_OK = .FALSE.
  ELSE; PRINT '(A)', 'PASSED'; END IF

  !=====================================================================
  ! Test SP 1: in-place R2C + C2R round-trip (N=30)
  !=====================================================================
  PRINT '(A)', '--- SP Test 1: in-place round-trip N=30 ---'
  DO I = 1, N1
    TS = 2.0_4*PIS*REAL(I-1,KIND=4)/REAL(N1,KIND=4)
    ORIGS1(I) = SIN(3.0_4*TS) + 0.5_4*COS(5.0_4*TS)
    RS1(I) = ORIGS1(I)
  END DO

  CALL sfftw_plan_dft_r2c_1d(PLAN_FWD,N1,RS1,CS1,FFTW_ESTIMATE)
  CALL sfftw_execute_dft_r2c(PLAN_FWD,RS1,CS1)
  CALL sfftw_destroy_plan(PLAN_FWD)

  CALL sfftw_plan_dft_c2r_1d(PLAN_BWD,N1,CS1,RS1,FFTW_ESTIMATE)
  CALL sfftw_execute_dft_c2r(PLAN_BWD,CS1,RS1)
  CALL sfftw_destroy_plan(PLAN_BWD)

  ERRS = 0.0_4
  DO I = 1, N1
    ERRS = MAX(ERRS,ABS(RS1(I)-ORIGS1(I)*REAL(N1,KIND=4)))
  END DO
  PRINT '(A,ES12.4)', 'Max error = ', REAL(ERRS,KIND=8)
  IF (ERRS > 1.0E-5_4) THEN; PRINT '(A)', 'FAILED'; ALL_OK = .FALSE.
  ELSE; PRINT '(A)', 'PASSED'; END IF

  !=====================================================================
  ! Test SP 2: out-of-place R2C + C2R round-trip (N=48)
  !=====================================================================
  PRINT '(A)', '--- SP Test 2: out-of-place round-trip N=48 ---'
  DO I = 1, N2
    TS = 2.0_4*PIS*REAL(I-1,KIND=4)/REAL(N2,KIND=4)
    ORIGS2(I) = COS(2.0_4*TS) - 0.3_4*SIN(7.0_4*TS)
    AS2(I) = ORIGS2(I)
  END DO

  CALL sfftw_plan_dft_r2c_1d(PLAN_FWD,N2,AS2,CS2,FFTW_ESTIMATE)
  CALL sfftw_execute_dft_r2c(PLAN_FWD,AS2,CS2)
  CALL sfftw_destroy_plan(PLAN_FWD)

  CALL sfftw_plan_dft_c2r_1d(PLAN_BWD,N2,CS2,BS2,FFTW_ESTIMATE)
  CALL sfftw_execute_dft_c2r(PLAN_BWD,CS2,BS2)
  CALL sfftw_destroy_plan(PLAN_BWD)

  ERRS = 0.0_4
  DO I = 1, N2
    ERRS = MAX(ERRS,ABS(BS2(I)-ORIGS2(I)*REAL(N2,KIND=4)))
  END DO
  PRINT '(A,ES12.4)', 'Max error = ', REAL(ERRS,KIND=8)
  IF (ERRS > 1.0E-5_4) THEN; PRINT '(A)', 'FAILED'; ALL_OK = .FALSE.
  ELSE; PRINT '(A)', 'PASSED'; END IF

  !=====================================================================
  ! Test SP 3: batch out-of-place (N=30, LOT=4)
  !=====================================================================
  PRINT '(A)', '--- SP Test 3: batch out-of-place N=30, LOT=4 ---'
  DO J = 1, LOT
    DO I = 1, N1
      TS = 2.0_4*PIS*REAL(I-1,KIND=4)/REAL(N1,KIND=4)
      ORIGS_BATCH(I,J) = SIN(REAL(J,KIND=4)*TS)
      ASBATCH(I,J) = ORIGS_BATCH(I,J)
    END DO
  END DO

  CALL sfftw_plan_many_dft_r2c(PLAN_FWD,1,[N1],LOT, &
       ASBATCH,[N1],1,N1, &
       CSBATCH,[N1/2+1],1,N1/2+1, &
       FFTW_ESTIMATE)
  CALL sfftw_execute_dft_r2c(PLAN_FWD,ASBATCH,CSBATCH)
  CALL sfftw_destroy_plan(PLAN_FWD)

  CALL sfftw_plan_many_dft_c2r(PLAN_BWD,1,[N1],LOT, &
       CSBATCH,[N1/2+1],1,N1/2+1, &
       ASBATCH,[N1],1,N1, &
       FFTW_ESTIMATE)
  CALL sfftw_execute_dft_c2r(PLAN_BWD,CSBATCH,ASBATCH)
  CALL sfftw_destroy_plan(PLAN_BWD)

  ERRS = 0.0_4
  DO J = 1, LOT
    DO I = 1, N1
      ERRS = MAX(ERRS,ABS(ASBATCH(I,J)-ORIGS_BATCH(I,J)*REAL(N1,KIND=4)))
    END DO
  END DO
  PRINT '(A,ES12.4)', 'Max error = ', REAL(ERRS,KIND=8)
  IF (ERRS > 1.0E-5_4) THEN; PRINT '(A)', 'FAILED'; ALL_OK = .FALSE.
  ELSE; PRINT '(A)', 'PASSED'; END IF

  !=====================================================================
  ! Test SP 4: batch in-place (N=30, LOT=4)
  !=====================================================================
  PRINT '(A)', '--- SP Test 4: batch in-place N=30, LOT=4 ---'
  DO J = 1, LOT
    DO I = 1, N1
      TS = 2.0_4*PIS*REAL(I-1,KIND=4)/REAL(N1,KIND=4)
      ORIGS_BATCH(I,J) = COS(REAL(J,KIND=4)*TS)
      ASBIP(I,J) = ORIGS_BATCH(I,J)
    END DO
  END DO

  CALL sfftw_plan_many_dft_r2c(PLAN_FWD,1,[N1],LOT, &
       ASBIP,[N1+2],1,N1+2, &
       CSBIP,[N1/2+1],1,N1/2+1, &
       FFTW_ESTIMATE)
  CALL sfftw_execute_dft_r2c(PLAN_FWD,ASBIP,CSBIP)
  CALL sfftw_destroy_plan(PLAN_FWD)

  CALL sfftw_plan_many_dft_c2r(PLAN_BWD,1,[N1],LOT, &
       CSBIP,[N1/2+1],1,N1/2+1, &
       RSBIP,[N1+2],1,N1+2, &
       FFTW_ESTIMATE)
  CALL sfftw_execute_dft_c2r(PLAN_BWD,CSBIP,RSBIP)
  CALL sfftw_destroy_plan(PLAN_BWD)

  ERRS = 0.0_4
  DO J = 1, LOT
    DO I = 1, N1
      ERRS = MAX(ERRS,ABS(RSBIP(I,J)-ORIGS_BATCH(I,J)*REAL(N1,KIND=4)))
    END DO
  END DO
  PRINT '(A,ES12.4)', 'Max error = ', REAL(ERRS,KIND=8)
  IF (ERRS > 1.0E-5_4) THEN; PRINT '(A)', 'FAILED'; ALL_OK = .FALSE.
  ELSE; PRINT '(A)', 'PASSED'; END IF

  !=====================================================================
  ! Test 9: memory allocation / free (smoke test)
  !=====================================================================
  PRINT '(A)', '--- Memory allocation smoke test ---'
  HDP = fftw_alloc_complex(100_C_SIZE_T)
  HSP = fftwf_alloc_complex(100_C_SIZE_T)
  IF (C_ASSOCIATED(HDP) .AND. C_ASSOCIATED(HSP)) THEN
    CALL fftw_free(HDP)
    CALL fftwf_free(HSP)
    PRINT '(A)', 'PASSED'
  ELSE
    PRINT '(A)', 'FAILED'
    ALL_OK = .FALSE.
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

END PROGRAM TEST_FFTW3
