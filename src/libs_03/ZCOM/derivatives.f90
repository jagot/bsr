!===================================================================
  SUBROUTINE DERIVATIV1 (N,H,FI,F1)
!===================================================================
! 1st order derivatives with the three-point formulas. 
! Linear extrapolations are made at the boundaries.  
! FI: input f(x); H: interval; F1: f'
!===================================================================

  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I
  REAL(8), INTENT (IN)  :: H
  REAL(8), INTENT (IN)  :: FI(N+1)
  REAL(8), INTENT (OUT) :: F1(N+1)
  REAL(8) :: H2

! f' from three-point formulas

  H2 = 2.d0*H
  DO I = 2, N
    F1(I) = (FI(I+1)-FI(I-1))/H2
  END DO

! Linear extrapolation for the boundary points

  F1(1) = 2.d0*F1(2)-F1(3)
  F1(N+1) = 2.d0*F1(N)-F1(N-1)

  END SUBROUTINE DERIVATIV1

!===================================================================
  SUBROUTINE DERIVATIV2 (N,H,FI,F2)
!===================================================================
! 2st order derivatives with the three-point formulas. 
! Linear extrapolations are made at the boundaries.  
! FI: input f(x); H: interval; F2: f''
!===================================================================

  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I
  REAL(8), INTENT (IN)  :: H
  REAL(8), INTENT (IN)  :: FI(N+1)
  REAL(8), INTENT (OUT) :: F2(N+1)
  REAL(8) :: HH

! f'' from three-point formulas

  HH = H*H
  DO I = 2, N
    F2(I) = (FI(I+1)-2.d0*FI(I)+FI(I-1))/HH
  END DO

! Linear extrapolation for the boundary points

  F2(1) = 2.d0*F2(2)-F2(3)
  F2(N+1) = 2.d0*F2(N)-F2(N-1)

  END SUBROUTINE DERIVATIV2

!===================================================================
  SUBROUTINE DERIVATIVE1_NONUNIFORM (N,X,F,F1)
!===================================================================
! 1st order derivatives with the three-point non-uniform formula. 
! Linear extrapolations are made at the boundaries.  
! X: input x; F: input f(x);  F1: f'
!===================================================================

  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I
  REAL(8), INTENT (IN)  :: X(N),F(N)
  REAL(8), INTENT (OUT) :: F1(N)
  REAL(8) :: h,hh, h1, hh1

! f' from three-point non-uniform formulas

  h1 = X(2)-X(1);  hh1=h1*h1
  DO i = 2, N-1
    h = X(i+1)-X(i); hh = h*h
    F1(i) = (hh1*F(i+1)+(hh-hh1)*F(i)-hh*F(i-1))/(h*h1*(h+h1))
    h1 = h; hh1 = hh
  END DO

! Linear extrapolation for the boundary points

  F1(1) = 2.d0*F1(2)-F1(3)
  F1(N) = 2.d0*F1(N-1)-F1(N-2)

  END SUBROUTINE DERIVATIVE1_NONUNIFORM

