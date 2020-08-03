!====================================================================
      Real(8) Function HWF(N,L,Z,R)
!====================================================================
!     returns the value of an unnormalized hydrogenic nl-function
!     with nuclear charge Z and radius R
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: N,L
      Real(8), intent(in) :: Z,R
      Integer(4) :: K, I
      Real(8) :: A,B,C,P,X

      K = N - L - 1
      if(K.lt.0) then
       write(*,'(a,a,i4,a,i4)') &
       '   HWF --> forbidden combination of n and l:', &
       '   N = ',N,'   L =',L
       Stop
      end if

      P = 1.d0
      A = 1.d0
      B = K
      C = N + L
      X = -2.d0*Z*R/N

      if(X.lt.-80.0) then          ! test if underflow may occur,
         HWF=0.0                   ! if so set hwf = 0
      else
         Do I=1,K
          P = 1.d0 + A/B*P/C*X
          A = A + 1.d0
          B = B - 1.d0
          C = C - 1.d0
         End do
         HWF = P*EXP(X/2.d0)*(-X)**(L+1)
      end if

      End Function HWF


!====================================================================
      Real(8) Function HNORM (N,L,Z)
!====================================================================
!     returns the value of the normalization constant for an (nl)
!     hydrogenic function with nuclear charge Z
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: N,L
      Real(8), intent(in) :: Z
      Integer :: M, I
      Real(8) :: A,B,D,T

      M = L + L + 1
      A = N + L
      B = M
      T = A
      D = B
      M = M - 1

      Do I=1,M
       A = A - 1.d0
       B = B - 1.d0
       T = T*A
       D = D*B
      End do

      HNORM = SQRT(Z*T)/(N*D)

      End Function HNORM


!======================================================================
      Subroutine HVA(n,l,z,ar1,ar2,ar3,am1,am2,am3)
!======================================================================
!     Avarage meaning of the r^1, r^2, r^3, r-1, r-2, r-3
!     for H-like functions
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,l
      Real(8), intent(in) :: z
      Real(8), intent(out) :: ar1,ar2,ar3,am1,am2,am3
      Integer :: nn, ll

      nn= n*n
      ll= l*(l+1)

      ar1 = 3*nn-ll
      ar1 = ar1/2.d0/z

      ar2 = 5*nn+1-3*ll
      ar2 = ar2*nn/2.d0/z**2

      ar3 = 35*nn*(nn-1) - 30*nn*(l+2)*(l-1) + 3*(l+2)*(l+1)*l*(l-1)
      ar3 = ar3*nn / 8.d0 /z**3

      am1 = z/nn

      am2 = z**2/nn/n/(l+0.5)

      am3 = 0.d0
      if(l.gt.0) am3 = z**3 / (nn*n*(l+1)*l*(l+0.5))

      End Subroutine HVA

