!=======================================================================
      Subroutine DCWF (n,kappa,Z,E,NTP,R,P,Q)
!=======================================================================
!   This subroutine computes the  Dirac-Coulomb  bound-state orbital
!   radial wavefunction.
!
!   Input:
!
!      n          The (usual) principal quantum number
!      kappa      The relativistic angular quantum number
!      Z          The effective nuclear charge
!      R(1:NTP)   Radial grid
!
!   Output:
!
!      E          The Dirac-Coulomb Eigenenergy (E-mc^2)
!
!      P          r times the large component wavefunction of
!
!      Q          r times the small component wavefunction of
!
!   Call(s) to: CGAMMA
!
!   Written by Farid A Parpia, at Oxford    Last Update: 14 Oct 1992
!
!=======================================================================
!   For a given value of n > 0, there are 2*n-1 eigenfunctions:
!   n  with  kappa = -1,-2,...,-n
!   n-1 with kappa = 1,2,...,n-1
!
!   k = iabs(kappa);  gamma = sqrt[k^2-(alfa*Z)^2]
!
!   N = sqrt[n^2 - 2(n-k)(k-gamma)]
!
!   E(n,k) = c^2 / sqrt[1 + (alfa*Z)^2 / (gamma+n-k)^2 ]
!
!   x = 2*Z/N*r
!
!   P_nk(r) = sqrt[1+E(n,k)/c^2]  N(n,k)  e^-x/2  x^gamma *
!             [(N-kappa) F2  -  (n-k) F1]
!
!   Q_nk(r) = sqrt[1-E(n,k)/c^2]  N(n,k)  e^-x/2  x^gamma *
!             [(N-kappa) F2  +  (n-k) F1]
!
!   N(n,k) = 1/[N*G(2*gamma+1)]  *
!            sqrt [ [Z*G(2*gamma+1+n-k)] / [2(n-k)!(N-kappa)] ]
!
!   G(x) - GAMMA function
!
!   F2 = F(-n+k,2*gamma+1,x);  F1 = F(-n+k+1,2*gamma+1,x)
!
!   F(a,b,x) = 1 + a/b x^1/1!  + [a(a+1)]/[b(b+1)] x^2/2! + ...
!=======================================================================
      Use zconst, only: alpha, C => c_au

      Implicit real(8) (A-H,O-Z)
      Real(8) :: R(*), P(*), Q(*)
      Real(8) :: T1(n),T2(n)

! ... Check the input arguments:

      IF(n.le.0)     Stop 'DCWF: Principal quantum number < 0'
      IF(kappa.eq.0) Stop 'DCWF: Kappa quantum number = 0'
      IF(kappa.eq.n) Stop 'DCWF: Kappa quantum number = n'
      IF(kappa.gt.n) Stop 'DCWF: Kappa quantum number > n'
      IF(Z.le.0.d0)  Stop 'DCWF: Nuclear charge is too small, <= 0'
      IF(Z.gt.C)     Stop 'DCWF: Nuclear charge exceeds limit, c'

! ... Now determine all the parameters:

      fn = DBLE (n)
      fkappa = DBLE (kappa)
      k = IABS (kappa);  fk = DBLE (k)
      nr = n-k;  fnr = DBLE (nr)
      Za = Z*alpha
      gamma = SQRT (fk*fk-Za*Za)
      gg = gamma + gamma + 1.d0
      BIGN = SQRT (fn*fn-2.d0*fnr*(fk-gamma))
      EPS = 1.d0 /SQRT(1.d0+(Za/(gamma+fnr))**2)

! ... EPS is the total energy divided by C*C

      E = (1.d0-EPS)*C*C         !  E => E-mc^2

! ... normalization constant N(n,k):

      NRFAC=1; Do I=1,NR; NRFAC = NRFAC*I; End do

      ARGI = 0.d0
      ARGR = gg+FNR;    CALL CGAMMA (ARGR,ARGI,RGAMM1,DUMMY)
      ARGR = gg;        CALL CGAMMA (ARGR,ARGI,RGAMM2,DUMMY)

      FAC = - SQRT (RGAMM1)/(RGAMM2*SQRT (DBLE(NRFAC))) &
            * SQRT (Z/(2.d0*BIGN*BIGN*(BIGN-FKAPPA)))

! ... Ensure that the slope of the large-component function is
! ... positive at the origin:

      IF (KAPPA .GT. 0) FAC = -FAC

      FG = FAC * SQRT(1.d0+EPS)
      FF = FAC * SQRT(1.d0-EPS)

! ...  Now set up the coefficients of the confluent hypergeometric
! ...  functions  F(-NR+1,2*GAMMA+1;RHO)  and  F(-NR,2*GAMMA+1;RHO)
! ...  in the workspace arrays  TA  and  TB, respectively

      nt1=0; if(nr.gt.0) nt1=nr-1;  nt2=nr

      FAC = 1.d0;  FACN = FAC
      A   = -fnr;  AN1  = A + 1.d0;  AN2 = A
      B   =   gg;  BN   = B

      K = 0
    2 K = K+1
      FDEN = 1.d0/(FACN*BN)
      if (K .LE. nt1)      T1(K) = AN1*FDEN
      if (K .LE. nt2) then
                           T2(K) = AN2*FDEN
         A = A + 1.d0;     AN1  = AN1*(A+1.d0); AN2 = AN2*A
         B = B + 1.d0;     BN   = BN*B
         FAC = FAC + 1.d0; FACN = FACN*FAC
         go to 2
      end if

! ...  Now tabulate the function over the entire grid

      FAC = (Z+Z)/BIGN
      a1 = FNR; a2 = BIGN-FKAPPA
      Do i = 1,NTP
       x = FAC*R(i); y = x
       F1 = 1.d0; F2 = 1.d0
       K = 0
    3  K = K+1
       if (K .LE. nt1) F1 = F1+T1(K)*y
       if (K .LE. nt2) then;  F2 = F2+T2(K)*y; y=y*x; go to 3;  end if

       F1 = a1*F1;  F2 = a2*F2
       OVLFAC = EXP(-0.5d0*x) * x**gamma

       P(I) = FG*OVLFAC*(F1-F2)
       Q(I) = FF*OVLFAC*(F1+F2)

      End do

      End Subroutine DCWF


!======================================================================
      Real(8) Function E_dcwf(n,k,Z)
!======================================================================
!     energy of the Dirac-Coulomb orbital nk:
!
!     E(n,k) = c^2 / sqrt[1 + (alfa*Z)^2 / (gamma+n-k)^2 ]  - c^2
!
!     gamma = sqrt[k^2-(alfa*Z)^2]
!----------------------------------------------------------------------
      Use zconst, only: alpha, C => c_au

      Implicit real(8) (A-H,O-Z)

      Za = Z*alpha

      gamma = SQRT (k*k-Za*Za)

      E = 1.d0 /SQRT(1.d0+(Za/(gamma+n-k))**2)

      E_dcwf = -(1.d0-E)*C*C

      End Function E_dcwf


!========================================================================
      Subroutine DCME(n,kappa,z,ar1,ar2,am1,am2,am3)
!========================================================================
! ... provides radial moments <r^k> for Dirac-Coulomb wave functions
! ... expressions are taken from DRAKE HANDBOOK,2006
!------------------------------------------------------------------------
!     let  x = 2Z * r, then
!
!     <x^2> = 2*N^2*[(5*N^2-2*kappa^2)*R^2 + (1-gamma^2) - 3*kappa*R]
!
!     <x>   = -kappa + (3*N^2-kappa^2)*R
!
!     <x-1> = [n*gamma + (k-gamma)*k] / [2*gamma*N^3]
!
!     <x-2> = [kappa^2 * R] / [2*gamma^2*N^3*(2*gamma-sgn(kappa))]
!
!     <x-3> = [N^2 + 2*gamma^2*kappa^2 - 3* N^2*kappa*R] /
!             [4*N^5*gamma*(gamma^2-1)*(4*gamma^2-1)]
!
!     where  R = sqrt[1 - Z^2/(N^2*c^2)]
!            N = sqrt[n^2-2*(n-k)(k-gamma)]
!            gamma = sqrt[kappa^2-(Z/c)^2]
!            k = abs(kappa)
!------------------------------------------------------------------------
      USE zconst, ONLY: c_au

      Implicit none
      Integer, intent(in) :: n, kappa
      Real(8), intent(in) :: z
      Real(8), intent(out) :: ar1,ar2,am1,am2,am3
      Integer :: k, kk
      Real(8) :: Za, gamma, BigN, R, RR, bb, gg, x,xx,xxx

      k = iabs(kappa); kk = k*k
      Za = (z/c_au)**2
      gamma = sqrt(kk-Za); gg = gamma*gamma;
      BigN = sqrt(n*n-2*(n-k)*(k-gamma)); bb=BigN*BigN
      R = sqrt(1-Za/bb); RR = R*R

      ar2 = 4*bb*( (5*bb-2*kk)*RR + (1-gg) - 3*kappa*R )
      ar1 = -kappa + (3*bb-kk)*R
      am1 = (n*gamma + (k-gamma)*k) / (2*gamma*bb*BigN)
      am2 = kk*R / (2*gg*bb*BigN*(2*gamma-sign(1,kappa)))
      am3 = (bb + 2*gg*kk - 3*bb*kappa*R) / &
            (4*bb*bb*BigN*gamma*(gg-1)*(4*gg-1))

      x = 2*Z !/BigN;
      xx=x*x; xxx=xx*x

      ar2 = ar2 / xx
      ar1 = ar1 / x
      am1 = am1 * x
      am2 = am2 * xx
      am3 = am3 * xxx

      End Subroutine DCME


!======================================================================
      Subroutine  Exp_dcwf(n,k,Z,ap,bp,cp)
!======================================================================
!     Expectation values for Dirac Coulomb w.f.
!     after notes of Mohr (???)
!----------------------------------------------------------------------
      Use zconst

      Implicit none
      Integer, intent(in) :: n,k
      Real(8), intent(in) :: Z
      Real(8) :: ap(-4:4),bp(-4:4),cp(-4:4)
      Real(8) :: la, nr, dn,dk, g, a, E, x,y, p
      Integer :: ip

      dn = DBLE(n); dk = DBLE(k)
      g = Z/c_au; la=sqrt(dk*dk-g*g); nr=dn-abs(dk)
      a = g/sqrt( (nr+la)**2+g*g )
      E = sqrt(one-a*a)

      ap(0)=one
      bp(0)=E
      cp(0)=dk * a**2 /g

      ap(-1) = a**3 / g**2 * (dk**2 /la + nr)
      bp(-1) = a*a / g
      cp(-1) = dk * a**3 / (g*la)

      ap(-2) = two*dk*a**3 / (g*la) * (two*dk*E-one) / (four*la*la-one)
      bp(-2) = two*a**3 / (g*la) * (two*la**2-dk*E) / (four*la*la-one)
      cp(-2) = two*a**3 / la * (two*dk*E-one) / (four*la*la-one)

      x = two*a**3 / la / (four*la*la-one)
      y = la*la-one

      ap(-3) = x * ( three*dk*E*(dk*E-one)/y - one)
      bp(-3) = x * E *( three*(la*la-dk*E)/y - one)
      cp(-3) = x /g * ( three*g*g*E*(dk*E-one)/y - a*a*dk)

      x = (four*la**2-three**2)
      y = E*bp(-3) - ap(-3)

      ap(-4) = two/three/x * &
               (four*g*dk*y + (-six*dk*E+four*g*g+nine)*cp(-3))
      bp(-4) = two/x * (two*g*y + (-three*E+two*dk)*cp(-3))
      ap(-4) = two/three/x * &
               ((four*dk**2-nine)*y + two*g*(two*dk-three*E)*cp(-3))

      ip=-3; p = DBLE(ip)
             x = (four*la*la-p*p)
             y = (E*bp(ip)-ap(ip))
      ap(ip-1) = two / p /x  * ( four*g*dk*y+ &
                (two*dk*p*E+four*g*g+p*p)*cp(ip) )
      bp(ip-1) = two / x * ( two*g*y + (two*dk+p*E)*cp(ip) )
      cp(ip-1) = two / p /x  * ( (four*dk*dk-p*p)*y + &
                 two*g*(two*dk+p*E)*cp(ip) )

      End Subroutine  Exp_dcwf

