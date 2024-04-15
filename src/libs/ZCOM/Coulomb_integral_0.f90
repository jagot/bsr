!======================================================================
      Real(8) Function IZ0_lamda (lamda, k1,l1, k2,l2)   result(F)
!======================================================================
!     Coulomb integrals,  Z = 0
!     (after A.Burgess and C. Whelam, CPC, 47, 295 (1987)
!-----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: lamda,l1,l2
      Real(8), intent(in) :: k1,k2

      Integer ::  L, m1,m2, ml, i,j
      Real(8) :: S, S1, S2, SL, pi, x, k
      Real(8), external :: E_fun, F0_lamda_diag
      Real(8) :: PL(0:2001), QL(0:2001) ! dimensions agree with DGELENI
      Real(8), intrinsic :: GAMMA

      F =  E_fun (l2,l1,lamda); if(F.eq.0.d0) Return

      if(abs(k1-k2)/(k1+k2).lt.1.d-6) then
        F =  F0_lamda_diag (lamda, k1,l1,l2); Return
      end if

      pi = ACOS(-1.d0)
      F = F * (-1)**((lamda+l2-l1)/2) * sqrt(pi/(k1*k2))  &
              / 2.d0**(lamda+1) / GAMMA(0.5d0 + lamda)

      M1 = (l1+l2-lamda)/2
      M2 = (l1+l2+lamda)/2

! ... Q-functions:

      ml = m2; if(ml.lt.3) ml=10

      x  = (k1**2 + k2**2) / (2.d0*k1*k2)

      Call DLEGENI(x,ml,PL,QL,i)

! ... main summation:

      SL = 0.d0
      Do L = M1,M2
       S = 0.d0
       Do i = 0, lamda; j = lamda-i
        S1 = E_fun(l1,L,j); if(s1.eq.0.d0) Cycle
        S2 = E_fun(l2,L,i); if(s2.eq.0.d0) Cycle
        S = S + k1**j * (-k2)**i /(s1*s2)
       End do
       if(S.eq.0.d0) Cycle
       SL = SL + (L+L+1)*QL(L)*S
      End do

      F = F * SL

      End Function IZ0_lamda


!======================================================================
      Real(8) Function E_fun(a,b,c)
!======================================================================
!     E =  (a+b+c+1)! (1/2(a+b-c))! (1/2(a-b+c))! (1/2(-a+b+c))!
!          / (1/2(a+b+c))!  / (a+b-c)!
!
!     a + b + c  - even   plus   triangle relation  {a,b,c}
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: a,b,c
      Real(8) :: X, Y
      Integer :: i, k, m, N1,N2, M1,M2,M3
      Integer, external :: ITRA

      E_fun = 0.d0;  m = a+b+c; if(mod(m,2).ne.0) Return
      if(ITRA(a,b,c).eq.0) Return

      m1 = (a+b-c)/2; m2 = (a-b+c)/2;  m3 = (b+c-a)/2
      n1 = (a+b-c);   n2 = (a+b+c)/2

      X = m + 1
      DO i = 2,m
       k = 1
       IF(M1.GE.I) K=K+1
       IF(M2.GE.I) K=K+1
       IF(M3.GE.I) K=K+1
       IF(N1.GE.I) K=K-1
       IF(N2.GE.I) K=K-1
       Y=I;  X = X * Y**K
      End do
      E_fun = X

      End Function E_fun


!======================================================================
      Real(8) Function F0_lamda_diag (lamda, k,l1,l2)
!======================================================================
!     F0_diag = pi k**(lamda-1) GAMMA(lamda) GAMMA[1/2(l1+l2-lamda+2)]/
!               2**(lamda+1) GAMMA[1/2(l1-l2+lamda+1)] GAMMA[1/2(l2-l1+lamda+1)]
!               GAMMA[1/2(l1+l2+lamda+2)]
!-----------------------------------------------------------------------

      Implicit none
      Integer, intent(in) :: lamda,l1,l2
      Real(8), intent(in) :: k

      Real(8) :: F, pi

      pi =  ACOS(-1.d0)

      F0_lamda_diag = pi *k**(lamda-1) /2.d0**(lamda+1) * GAMMA(1.d0*lamda) &
                      / GAMMA((l1-l2+lamda+1.d0)/2) * GAMMA((l1+l2-lamda+2.d0)/2) &
                      / GAMMA((l2-l1+lamda+1.d0)/2) / GAMMA((l1+l2+lamda+2.d0)/2)

      End Function F0_lamda_diag


!=============================================================================
      SUBROUTINE DLEGENI(Z,NMAX,PL,QL,NUEVO)
!=============================================================================
!     LEGENDRE FUNCTIONS OF FIRST AND SECOND KIND OF ARGUMENT GREATER THEN ONE
!     (after  A.Gil and J.Segura, CPC, 105, 273, 1997)
!-----------------------------------------------------------------------------
!  INPUT :
!
!    Z        ARGUMENT OF THE FUNCTIONS
!
!    NMAX     MAXIMUM ORDER OF THE FUNCTIONS :
!             WE SHALL GET  FUNCTIONS OF ALL THE ORDERS BELOW
!             MIN(NMAX). NUEVO IS DEFINED BELOW .
!
!  OUTPUT :
!
!    PL(L+1)
!             WE SHALL KEEP THESE VALUES IN AN ARRAY
!    QL(L+1)
!             WE SHALL KEEP THESE VALUES IN AN ARRAY
!    NUEVO    MAXIMUM ORDER OF  FUNCTIONS CALCULATED WHEN
!             PL (NMAX+1,Z) IS LARGER THAN 10**NPRE
!
!---------------------------------------------------------------------------

        IMPLICIT REAL*8 (A-H,O-Z)

        Real(8) :: PL(0:2001),QL(0:2001)
        Real(8), parameter :: PI=3.14159265358979323D0
        Real(8), parameter :: EPS=1.D-16
        Real(8), parameter :: TINY=1.D-180
        Integer, parameter :: NPRE=300

!   EPS: REQUIRED ACCURACY FOR THE CONTINUED FRACTION (MODIFIED LENTZ)
!   TINY: SMALL PARAMETER TO PREVENT OVERFLOWS

        IF (Z.LE.1.D0) THEN
          WRITE(6,*)'IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1'
          STOP
        END IF

!   WE USE THE CODE IF NMAX IS GREATER THAN OR EQUAL TO 2

        NMAXP=NMAX
        IF(NMAX.LT.2) NMAX=2

!   LIMITS TO THE VALUE OF NMAX:
!   THE MAXIMUM VALUE OF NMAX IS SUCH THAT
!   PL(NMAX+1) IS LESS THAN 10**NPRE

         ALFA=DLOG(Z+DSQRT(Z*Z-1.D0))

         CANT=0.5*( DLOG(2.D0*PI)+DLOG(DSINH(ALFA))-ALFA)+NPRE*DLOG(10.D0)

         DO WHILE ((ALFA*NMAX-0.5D0*DLOG(DFLOAT(NMAX))).GE.(CANT))
           NMAX=NMAX-3
         END DO

         NUEVO=NMAX

!   WE STORAGE THE VALUES OF PL(0) AND PL(1)

         PL(0)=1.D0
         PL(1)=Z

!   WE USE THE RECURRENCE RELATIONS

        DO NP=1,NMAX
          PL(NP+1)=((2.D0*NP+1.D0)*Z*PL(NP)-NP*PL(NP-1))/(NP+1.D0)
        END DO

!   WE EVALUATE THE CONTINUED FRACTION USING
!   LENTZ ALGORITHM

        N=NMAX
        M=0
        B=(2.D0+1.D0/(N+1.D0) )*Z
        A=1.D0
        FC=TINY
        C0=FC
        D0=0.D0
10      D0=B+A*D0
        IF (D0.EQ.0.D0) D0=TINY
        C0=B+A/C0
        IF (C0.EQ.0.D0) C0=TINY
        D0=1.D0/D0
        DELTA=C0*D0
        FC=FC*DELTA
        M=M+1
        A=-(1.D0+1.D0/DFLOAT(N+M))
        B=(2.D0+1.D0/(N+M+1.D0))*Z
        IF (ABS(DELTA-1.D0).GT.EPS) GOTO 10

!     WE EVALUATE QL(NMAX+1) AND QL(NMAX) USING :
!     THE WRONSKIAN : W{PL(NMAX),QL(NMAX)} =1./(1.-Z**2)
!     THE KNOWN VALUES OF PL(NMAX+1) AND PL(NMAX)
!     THE VALUE OF H = QL(NMAX+1)/QL(NMAX)

        QL(NMAX)=1.D0/(PL(NMAX+1)-FC*PL(NMAX))/(NMAX+1.D0)
        QL(NMAX+1)=QL(NMAX)*FC

        DO I=1,NMAX
          NP=NMAX+1-I
          N=NP-1
          QL(N)=((NP+NP+1.D0)*Z*QL(NP)-(NP+1.D0)*QL(NP+1))/NP
        END DO

        NMAX=NMAXP

        END SUBROUTINE DLEGENI
