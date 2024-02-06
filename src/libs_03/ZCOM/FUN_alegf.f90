!======================================================================
      Real(8) Function ALEGF (L,M,X,norm)
!======================================================================
! ... ASSOCIATED LEGENDRE POLYNOMIALS:
! ... NORM=0: UNNORMALISED LEGENDRE FUNCTIONS P_LM(x) (I.E. ORIGINAL)
! ... NORM=1: NORMALISED: * SQRT[ (2L+1)/2 * (L-M)!/(L+M)! ]
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Real(8) :: zero = 0.d0; one = 1.d0; two = 2.d0

      ALEGF=zero; MA=IABS(M); if(L.lt.MA) Return

      ALEGF=one
      FACT=one

      if(M.ne.0) then
       Do I=L-MA+1,L+MA; FACT=FACT*I; End do; FACT=one/FACT
      end if

      if(M.lt.0) ALEGF=FACT
      if(norm.eq.1.and.m.ge.0) ALEGF = ALEGF * sqrt(FACT*(L+L+1)/2)
      if(norm.eq.1.and.m.lt.0) ALEGF = ALEGF * sqrt((L+L+1)/(2*FACT))
      if(L.eq.0) Return


      if(abs(x).eq.one) then
       if(m.ne.0) ALEGF=zero; if(x.lt.zero) ALEGF=ALEGF*(-1)**L
       Return
      end if


      FACT=one; Do I=1,L+L-1,2; FACT=FACT*I; End do

      if(MA.eq.L) then
       ALEGF=ALEGF*FACT*(one-x*x)**(0.5d0*MA)
       Return
      end if
      if(MA.eq.L-1) then
       ALEGF=ALEGF*FACT*(one-x*x)**(0.5d0*MA)*x
       Return
      end if

! ... for m <= 1, use the recurrence relations in respect to l:
! ... P(l+1) = [(2l+1)P(l)-(l+m)p(l-1)] / (l-m+1)

      if(MA.le.1) then
       if(MA.eq.0) then
        P0=one;  P1=x
       else
        P0=zero; P1=DSQRT(one-x*x)
       end if
       DO N=1,L-1
        PN=((N+N+1)*x*P1-(N+MA)*P0)/(N-MA+1); P0=P1; P1=PN
       End do
        ALEGF=PN*ALEGF; Return
      end if

! ... use the recurrence relations in respect to m:
! ... P(m+2) = 2(m+1)x/sqrt(1-x^2)P(m+1)-(l-m)(l+m+1)P(m)

      z = x/sqrt(one-x*x)
      P2 = FACT*(one-x*x)**(0.5d0*L)
      P1 = P2 * z; z = z * two
      Do MM = L-2,MA,-1
       PM=((MM+1)*z*P1-P2)/((L-MM)*(L+MM+1)); P2=P1; P1=PM
      End do

      ALEGF=PM*ALEGF; Return

      End Function ALEGF



!=====================================================================
      Real(8) Function ALEGFM (L,M,X,norm)
!=====================================================================
!     sign correction for origional function ALEGF
!---------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      ALEGFM = ALEGF (L,M,X,norm)
      if(m.gt.0) ALEGFM = ALEGFM * (-1)**M
      End Function ALEGFM

