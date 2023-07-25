!----------------------------------------------------------------------
      SUBROUTINE P_derive(i,j)
!----------------------------------------------------------------------
!
!     Define the derivative of the P_i: 
!
!     (d - 1/r) P_i   -->  (r^1/2)(d - 1/2r)) F_i  --> P_j
!
!----------------------------------------------------------------------

      USE RADIAL, L => lro, N => nro, KS => kro, EL => ero, MX => mro

      IMPLICIT NONE
      Integer(4), INTENT(in) :: i,j      
      Integer(4) :: m
      Real(8), Allocatable, Dimension(:) :: A,B,C,D
      Real(8) :: S

      if(.not.allocated(A)) Allocate(A(nr),B(nr),C(nr),D(nr))

      if(nrf.gt.mrf-2) Call Alloc_radial(mrf+jrf)

      A = P(:,i);   Call SPLIN3 (nr, R, A,B,C,D)

      P(:,j) =  B - D5*A/R

      P(mx(i)+1:NR,j) = D0

      N(j) = N(i)
      L(j) = L(i)
      KS(j) = KS(i)
      EL(j) = EL(i)
      MX(j) = MX(i)

      if(L(j).gt.0) then
       AZ(j) = L(i)*AZ(i)
       mexp(j) = L(j)
       aexp(j) = -D1/L(j)
       bexp(j) = bexp(i)*(L(j)+2)/L(j)
      else
       AZ(j) = -Z*AZ(i)
       mexp(j) = 1
       aexp(j) = - D2*bexp(i)
       bexp(j) = 0.d0
       bexp(j) = P(4,j)/(AZ(j)*R2(4))
       S = Z*R(4)
       bexp(j) = (bexp(j) - D1 - S*aexp(j))/S**2
      end if

      Do m = 1,3
       S = Z*R(m)
       P(m,j) = AZ(j)*r(m)**mexp(j)*(D1 + S*(aexp(j) + S*bexp(j)))
       P(m,j) = P(m,j)/R2(m)
      End do

      END SUBROUTINE P_derive



!======================================================================
      SUBROUTINE SPLIN3 (N, X, Y, B, C, D)
!======================================================================
!
!     THE COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N  ARE COMPUTED
!     FOR A CUBIC INTERPOLATING SPLINE
!
!     S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
!
!     FOR  X(I) .LE. X .LE. X(I+1)
!
!   INPUT ...
!
!     N = THE NUMBER OF DATA POINTS OR KNOTS (N.GE.2)
!     X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
!     Y = THE ORDINATES OF THE KNOTS
!
!   OUTPUT ...
!
!     B, C, D  = ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.
!
!     USING  P  TO DENOTE DIFFERENTIATION:
!     Y(I) = S(X(I))
!     B(I) = SP(X(I))
!     C(I) = SPP(X(I))/2
!     D(I) = SPPP(X(I))/6  (DERIVATIVE FROM THE RIGHT)
!
!     THE ACCOMPANYING FUNCTION  'SEVAL'  CAN BE USED
!     TO EVALUATE THE SPLINE.
!
!----------------------------------------------------------------------
      
      IMPLICIT REAL(8) (A-H,O-Z)
      
      DIMENSION X(*),Y(*),B(*),C(*),D(*)
 
      NM1 = N-1
      IF ( N .LT. 2 ) RETURN
      IF ( N .LT. 3 ) GO TO 50
 
! ... SET UP TRIDIAGONAL SYSTEM

! ... B = DIAGONAL, D = OFFDIAGONAL, C = RIGHT HAND SIDE.
 
      D(1) = X(2) - X(1)
      C(2) = (Y(2) - Y(1))/D(1)
      Do I = 2, NM1
         D(I) = X(I+1) - X(I)
         B(I) = 2.*(D(I-1) + D(I))
         C(I+1) = (Y(I+1) - Y(I))/D(I)
         C(I) = C(I+1) - C(I)
      End do
 
! ... END CONDITIONS. THIRD DERIVATIVES AT  X(1)  AND  X(N)
! ... OBTAINED FROM DIVIDED DIFFERENCES
 
      B(1) = -D(1)
      B(N) = -D(N-1)
      C(1) = 0.
      C(N) = 0.
      IF ( N .EQ. 3 ) GO TO 15
      C(1) = C(3)/(X(4)-X(2)) - C(2)/(X(3)-X(1))
      C(N) = C(N-1)/(X(N)-X(N-2)) - C(N-2)/(X(N-1)-X(N-3))
      C(1) = C(1)*D(1)**2/(X(4)-X(1))
      C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))
 
!     FORWARD ELIMINATION
 
   15 Do I = 2, N
         T = D(I-1)/B(I-1)
         B(I) = B(I) - T*D(I-1)
         C(I) = C(I) - T*C(I-1)
      End do
 
! ... BACK SUBSTITUTION
 
      C(N) = C(N)/B(N)
      Do IB = 1, NM1
         I = N-IB
         C(I) = (C(I) - D(I)*C(I+1))/B(I)
      End do
 
! ... C(I) IS NOW THE SIGMA(I) OF THE TEXT
 
! ... COMPUTE POLYNOMIAL COEFFICIENTS
 
      B(N) = (Y(N) - Y(NM1))/D(NM1) + D(NM1)*(C(NM1) + 2.*C(N))
      Do I = 1, NM1
         B(I) = (Y(I+1) - Y(I))/D(I) - D(I)*(C(I+1) + 2.*C(I))
         D(I) = (C(I+1) - C(I))/D(I)
         C(I) = 3.*C(I)
      End do
      C(N) = 3.*C(N)
      D(N) = D(N-1)
      RETURN
!-----------------------------------------------------------------------
!                                       liniar interpolation for n = 2 :
   50 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.
      D(1) = 0.
      B(2) = B(1)
      C(2) = 0.
      D(2) = 0.

      RETURN
      END

