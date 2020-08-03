!======================================================================
      Subroutine SPLIN3 (N, X, Y, B, C, D)
!======================================================================
!     DEFINES COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N  
!     FOR A CUBIC INTERPOLATING SPLINE
!
!     S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
!
!     FOR  X(I) .LE. X .LE. X(I+1)
!
!     INPUT:
!
!     N = THE NUMBER OF DATA POINTS OR KNOTS (N=>2)
!     X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
!     Y = THE ORDINATES OF THE KNOTS
!
!     OUTPUT:
!
!     B, C, D  -  ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.
!
!     Y(I) = S(X(I))
!     B(I) = S'(X(I))
!     C(I) = S''(X(I))/2
!     D(I) = S'''(X(I))/6  (DERIVATIVE FROM THE RIGHT)
!
!     ACCOMPANYING FUNCTION  'SEVAL'  CAN BE USED TO EVALUATE THE SPLINE.
!----------------------------------------------------------------------
      Implicit none
      Integer, Intent(in ) :: N
      Real(8), Intent(in ) :: X(n),Y(n)
      Real(8), Intent(out) :: B(n),C(n),D(n)
      Integer :: I,J,M
      Real(8) :: T
 
      IF ( N .LT. 2 ) RETURN

      IF ( N .LT. 3 ) THEN ! liniar interpolation for n = 2 :

       B(1) = (Y(2)-Y(1))/(X(2)-X(1)); C(1)=0.d0; D(1)=0.d0
       B(2) = B(1);  C(2)=0.d0; D(2)=0.d0
       RETURN

      END IF

      M = N-1
 
! ... SET UP TRIDIAGONAL SYSTEM

! ... B = DIAGONAL, D = OFFDIAGONAL, C = RIGHT HAND SIDE.
 
      D(1) = X(2) - X(1)
      C(2) = (Y(2) - Y(1))/D(1)
      Do I = 2, M
         D(I) = X(I+1) - X(I)
         B(I) = 2*(D(I-1) + D(I))
         C(I+1) = (Y(I+1) - Y(I))/D(I)
         C(I) = C(I+1) - C(I)
      End do
 
! ... END CONDITIONS. THIRD DERIVATIVES AT  X(1)  AND  X(N)
! ... OBTAINED FROM DIVIDED DIFFERENCES
 
      B(1) = -D(1)
      B(N) = -D(N-1)
      C(1) = 0.
      C(N) = 0.
      IF ( N .GT. 3 ) THEN
       C(1) = C(3)/(X(4)-X(2)) - C(2)/(X(3)-X(1))
       C(N) = C(N-1)/(X(N)-X(N-2)) - C(N-2)/(X(N-1)-X(N-3))
       C(1) = C(1)*D(1)**2/(X(4)-X(1))
       C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))
      END IF
 
!     FORWARD ELIMINATION
 
      Do I = 2, N
         T = D(I-1)/B(I-1)
         B(I) = B(I) - T*D(I-1)
         C(I) = C(I) - T*C(I-1)
      End do
 
! ... BACK SUBSTITUTION
 
      C(N) = C(N)/B(N)
      Do J = 1, M
       I=N-J; C(I) = (C(I) - D(I)*C(I+1))/B(I)
      End do
 
! ... C(I) IS NOW THE SIGMA(I) OF THE TEXT
 
! ... COMPUTE POLYNOMIAL COEFFICIENTS
 
      B(N) = (Y(N) - Y(M))/D(M) + D(M)*(C(M) + 2.*C(N))
      Do I = 1, M
         B(I) = (Y(I+1) - Y(I))/D(I) - D(I)*(C(I+1) + 2.*C(I))
         D(I) = (C(I+1) - C(I))/D(I)
         C(I) = 3.*C(I)
      End do
      C(N) = 3.*C(N);  D(N) = D(N-1)

      End Subroutine SPLIN3


!=======================================================================
      Real(8) Function SEVAL(N,U,X,Y,B,C,D)
!=======================================================================
!     Evaluates the value of cubic spline for abscisa U
!-----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: N
      Real(8), intent(in) :: U
      Real(8), intent(in) :: X(N),Y(N),B(N),C(N),D(N)
      Integer :: mflag, I
      Real(8) :: DX

      Call INTERV (X,N,U,I,mflag)
      
      DX=U-X(I); SEVAL=Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))

      End function SEVAL

    