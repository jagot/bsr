!======================================================================
      Real(8) Function Ikl (ekk1,ll1,ekk2,ll2,lamda,acc1,deltak1,istep)
!======================================================================
!     Coulomb integrals,  Z > 0
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Integer, intent(in) :: ll1,ll2,lamda,istep
      Real(8), intent(in) :: ekk1,ekk2,acc1
      Integer :: l1,l2, jstep=1000
      Real(8) :: acc = 1.d-12, deltak1, deltak=1.e-5, k1,k2
      Real(8), external :: Ckl_factor

      Ikl = 0.d0       ! sign of failor

! ... check input parameters:

      if(abs(ll1-ll2).gt.lamda) Return
      if(abs(ll1+ll2).lt.lamda) Return
      if(ekk1+ekk2.le.0.d0) Return

! ... check the tolerence:

      if(acc1.gt.0.d0) acc = acc1
      if(deltak1.gt.0.d0) deltak = deltak1
      if(istep.gt.0) jstep=istep

! ... ek1 = ek2

      if(abs(ekk1-ekk2).le.deltak) then
        ek = ekk1; l = max(ll1,ll2)            !  dipole case only ???
        IKL = sqrt(ek)/sqrt(1.d0 + l*l*ek)/2.d0
        Return
      end if

! ... ek1 \= ek2

      if(ekk1.ge.ekk2) then
        l1=ll1; l2=ll2; ek1=ekk1; ek2=ekk2
      else
        l1=ll2; l2=ll1; ek1=ekk2; ek2=ekk1
      end if

      k1 = sqrt(ek1); k2 = sqrt(ek2)
      pi = acos(-1.d0)

!debug output:


      Call JLM(l1,l2-lamda,ek1,b0,b1)

      a0 = 1.d0
      SS = a0*b0
      S = S + ss

      a1 = -a0 /(l2+1)
      ss = a1*b1
      S = S + ss

      m  = l2 - lamda  + 1

      Do j = 2,jstep
       a2 = -(2.d0*a1 + ek2*a0) / (j*(2*l2+1+j))
       m = m + 1
       b2 = (-2.d0*b1 - (m+l1)*(m-l1-1)*b0) / ek1
       ss = a2*b2
       S = S + ss

       if(abs(ss).lt.acc*abs(S).and.j.gt.10) Exit

       a0 = a1; a1 = a2
       b0 = b1; b1 = b2

       b = max(abs(b0),abs(b1))
       if(b.gt.1.d10) then
        b0=b0/b; b1=b1/b;  a0 = a0*b;  a1 = a1*b
       end if

      End do

      Ikl = S * Ckl_factor(l2,ek2,1)

      End Function Ikl


!======================================================================
      Real(8) Function IKL0 (l,ekk1,ekk2,bcc,istep)
!======================================================================
!     Coulomb integrals, lamda=0
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Integer, intent(in) :: l,istep
      Real(8), intent(in) :: ekk1,ekk2, bcc
      Integer :: l1,l2, jmax=1000000  ! ???
      Real(8) :: eps = 1.d-12, k1,k2
      Real(8), external :: Ckl_factor, JLLM

! ... equal energies:

      IKL0 = 0.d0
      if(ekk1.eq.ekk2) Return

! ... accuracy limits:

      if(bacc.gt.0.d0) acc = bcc
      if(istep.gt.0) jstep=istep

! ... choose lower energy:

      if(ekk1.gt.ekk2) then
        ek1=ekk1; ek2=ekk2
      else
        ek1=ekk2; ek2=ekk1
      end if

      k1 = sqrt(ek1); k2 = sqrt(ek2)
      pi = acos(-1.d0)
      l1=l; l2=l

      Call JLLP (l,l,ek1,b0,b1)

      a0 = 1.d0
      SS = a0*b0
      S = S + ss

      a1 = -a0 /(l2+1)
      ss = a1*b1
      S = S + ss;   ss1=ss

      m  = l + 1

      Do j = 2,jmax
       a2 = -(2.d0*a1 + ek2*a0) / (j*(2*l2+1+j))
       m = m + 1
       b2 = (-2.d0*b1 - (m+l1)*(m-l1-1)*b0) / ek1
       ss = a2*b2
       S = S + ss

       if(abs(ss)+abs(ss1).lt.eps*abs(S)) Exit
       ss1 = ss

       a0 = a1; a1 = a2
       b0 = b1; b1 = b2

       b = max(abs(b0),abs(b1))
       if(b.gt.1.d10) then
        b0=b0/b; b1=b1/b;  a0 = a0*b;  a1 = a1*b
       end if

      End do

      Ikl0 = S * Ckl_factor(l2,ek2,1)

      End Function Ikl0



!======================================================================
      Real(8) Function Ckl_factor(l,ek,met)
!======================================================================
!     l   -  orbital momentum
!     ek  -  k*k, scaled energy in Ry
!            where  k  -  scaled liniar momemut,  k/z
!
!     C(k,l) =  (2pi)^1/2 * 2^l /  (2l+1)!
!             * [ Product(i=0,l) {i^2 k^2 +1)] ^ 1/2
!             / [1-exp(-2pi/k)] ^ 1/2
!
!     this fator is used fot the Coulomb-function normalization
!     with  F ~  1/k^1/2  sin(kr + ...)
!----------------------------------------------------------------------
      Implicit none

      Integer, intent(in) :: l, met
      Real(8), intent(in) :: ek

      Integer :: i
      Real(8) :: pi,k,S

      pi = acos(-1.d0)
      k = sqrt(ek)

      S = 2.d0*pi
      Do i = 1,l
       S = S * (1.d0 + i*i*ek)
      End do

      Ckl_factor = sqrt ( S/(1.d0 - exp(-2.d0*pi/k)) ) * 2.d0**l

      if(met.eq.-1) Return

      Do i = 2,l+l+1
        Ckl_factor = Ckl_factor  / i
      End do

      End Function Ckl_factor


!======================================================================
      Subroutine JLM (l,m,ek,J0,J1)
!======================================================================
!     Calculate the J(l,m) [J0] and J(l,m+1) [J1]  coefficients
!----------------------------------------------------------------------

      Implicit none
      Integer, intent(in) :: l,m
      Real(8), intent(in) :: ek
      Real(8), intent(out) :: J0,J1

      Integer :: i
      Real(8) ::  C, k, pi = 0.d0
      Real(8), external :: Ckl_factor, JLLM

      pi = acos(-1.d0)
      k  = sqrt(ek)

      if(m.ge.l) then

       J0 =  Ckl_factor(l,ek,-1) / ek**(l+1) * exp(-pi/k)  ! J_l^l
       J1 =  -2.d0 * J0 / ek                               ! J_l^(l+1)
       if(l.eq.m) Return
       Do i = l+1,m
        C = (-2.0*J1 - (i+l)*(i-l-1)*J0) / ek
        J0 = J1
        J1 = C
       End do
       Return

      else

       C = JLLM (l,ek,J1,J0)
       if(m.eq.l-1) Return
       Do i = l-2,m,-1                                 ! ???
        C = (-2.0*J0 - ek*J1) / (i+l+2)/(i-l+1)
        J1 = J0
        J0 = C
       End do

      end if

      End Subroutine JLM


!======================================================================
      Real(8) recursive Function JLLM (l,ek,JLL,JLM)  Result(S)
!======================================================================
!     Calculate the J(l,l) [JLL] and J(l,l-1) [JLM]  coefficients
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l
      Real(8) :: JLL, JLM, ek, pi,  C, ep
      Integer :: k

      pi = acos(-1.d0)

      if(l.eq.0) then

       ep  = exp(-pi/sqrt(ek))
       S   = sqrt( 2.d0 * pi / ( 1.d0 - ep*ep) )
       JLL = S * ep / ek
       JLM = S * ( 1.d0 - ep ) / 2.d0

      else

       S = JLLM(l-1,ek,JLL,JLM)
       C = sqrt(1.d0+l*l*ek)
       JLM = (l*(l+l-1)*JLM - JLL) / C
       JLL = JLL * 2.d0 * C / ek

      end if

      End Function JLLM



!======================================================================
      Subroutine JLLP (l,m,ek,JLL,JLP)
!======================================================================
      Implicit none
      Integer, intent(in) :: l,m
      Real(8), intent(in) :: ek
      Integer :: i
      Real(8) :: JLL,JLP, pi, C, k
      Real(8), external :: Ckl_factor

      pi = acos(-1.d0)
      k  = sqrt(ek)

      JLL =  Ckl_factor(l,ek,-1) / ek**(l+1) * exp(-pi/k)
      JLP =  -2.d0 * JLL / ek

      if(l.eq.m) Return

      Do i = l+1,m
       C = (-2.0_16*JLP - (i+l)*(i-l-1)*JLL) / ek
       JLL = JLP
       JLP = C
      End do

      End Subroutine JLLP











