PROGRAM DUMMY
  IMPLICIT NONE

  INTEGER, PARAMETER :: M = 2, N = 2, K = 3
  REAL(8) :: A(M,K), B(K,N), C(M,N)
  REAL(8) :: ALPHA, BETA
  INTEGER :: I, J
  EXTERNAL :: DGEMM

  ! A = [1 2 3; 4 5 6]
  A(1,1) = 1.0D0; A(1,2) = 2.0D0; A(1,3) = 3.0D0
  A(2,1) = 4.0D0; A(2,2) = 5.0D0; A(2,3) = 6.0D0

  ! B = [1 2; 3 4; 5 6]
  B(1,1) = 1.0D0; B(1,2) = 2.0D0
  B(2,1) = 3.0D0; B(2,2) = 4.0D0
  B(3,1) = 5.0D0; B(3,2) = 6.0D0

  ! C = A * B
  ! Expected: C = [22 28; 49 64]

  ALPHA = 1.0D0
  BETA  = 0.0D0
  C = 0.0D0

  CALL DGEMM('N', 'N', M, N, K, ALPHA, A, M, B, K, BETA, C, M)

  PRINT *, 'Result of DGEMM:'
  DO I = 1, M
    PRINT '(2F8.2)', (C(I,J), J = 1, N)
  END DO

  ! Check correctness
  IF (ABS(C(1,1)-22.0D0) < 1.0D-12 .AND. &
      ABS(C(1,2)-28.0D0) < 1.0D-12 .AND. &
      ABS(C(2,1)-49.0D0) < 1.0D-12 .AND. &
      ABS(C(2,2)-64.0D0) < 1.0D-12) THEN
    PRINT *, 'DGEMM test PASSED'
  ELSE
    PRINT *, 'DGEMM test FAILED'
    STOP 1
  END IF

END PROGRAM DUMMY
