!========================================================================
      SUBROUTINE DCME(n,kappa,z,ar1,ar2,am1,am2,am3)
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

      End SUBROUTINE DCME


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
