!====================================================================
      Module DBS_nuclear
!====================================================================
!     nuclear parameters
!--------------------------------------------------------------------
      Use zconst

      Implicit none

      Real(8) :: atomic_number  = zero
      Real(8) :: atomic_weight  = zero

      Character(20) :: nuclear = 'Fermi'

!     nuclear = point   - point nuclear
!             = uniform - uniform destribution
!             = Fermi   - Fermi destribution
!
!     uniform distribution:
!     V(r) = -3/2 Z/R [1-r^2/(3R^2],      r <  R
!          = -Z/r                         r >  R

      Real(8) :: r_uniform = zero
      Real(8) :: ro_uniform = zero

! ... Fermi distribution:   rho(r) = rho'/[1+exp((r-c)/a)]

      Real(8) :: a_fermi = zero
      Real(8) :: c_fermi = zero
      Real(8) :: t_fermi = 2.30   ! in fermi
      Real(8) :: rrms    = zero
      Real(8) :: ro_fermi = zero

! ... nuclear potential in different basis:

      Real(8), allocatable :: VR_nucl(:,:)
      Real(8), allocatable :: VB_nucl(:,:)
      Real(8), allocatable :: VF_nucl(:,:)
      Real(8), allocatable :: ZR_nucl(:,:)

! ... Nuclear spin (I) (in units of h/2*pi):

      Real(8) :: I_nuclear = zero

! ... Nuclear dipole moment (in nuclear magnetons):

      Real(8) :: D_nuclear = zero

! ... Nuclear quadrupole moment (in barns):

      Real(8) :: Q_nuclear = zero

      End Module DBS_nuclear


!======================================================================
      Subroutine Read_nuclear(nu,z,atw)
!======================================================================
!     Set up the nuclear parameters in module DBS_nuclear
!     z,awt - data from the calling program
!     if they are zero - try to get them from file unit "nu"
!----------------------------------------------------------------------
      Use DBS_nuclear

      Implicit none
      Integer :: nu, an
      Real(8) :: apar,cpar, z,atw, rms, A
      Character(200) :: atom, core, conf

! ... read parameters from file:

      if(nu.gt.0) then

       Call Read_rpar(nu,'atomic_number',atomic_number)
       Call Read_rpar(nu,'atomic_weight',atomic_weight)
       Call Read_apar(nu,'nuclear',nuclear)
       Call Read_rpar(nu,'r_uniform',r_uniform)
       Call Read_rpar(nu,'a_fermi',a_fermi)
       Call Read_rpar(nu,'c_fermi',c_fermi)
       Call Read_rpar(nu,'t_fermi',t_fermi)
       Call Read_rpar(nu,'rrms',rrms)
       if(z.le.0.d0) z = atomic_number
       atw = atomic_weight

      end if

      if(z.le.0.d0) Stop 'Stop in Read_nuclear:  Z <= 0'
      an = NINT(z)
      if(an.lt.1.or.an.gt.104) Stop 'Stop in Read_nuclear:  Z out of range'
      Call Def_atom(an,atom,A,rms,core,conf)

      if(atw.le.0.d0) atw = A
      if(rrms.eq.0.d0) rrms = rms
      if(rrms.lt.0.d0) then
       if(an.le.90)  then
         rrms=(0.836d0*atw**(1.d0/3.d0)+0.570d0)
       else
         rrms=(0.77d0*atw**(1.d0/3.d0)+0.980d0)
       end if
      end if

      if(nuclear.eq.'Fermi'.and.rrms.lt.2.d0) then
       rrms = 2.d0
       write(*,'(70("_"))')
       write(*,*) 'rrms is changed to 2.0 because GETCPR routine from GRASP'
       write(*,*) 'is nor working for small values of rrms'
       write(*,'(70("_")/)')
      end if

      atomic_number = z
      atomic_weight = atw

! ... guess for nucler radius:

      if(nuclear.eq.'uniform'.and.r_uniform.eq.0.d0) then
       r_uniform = rrms * fermi_in_cm / bohr_radius_in_cm
      end if

! ... check the Fermi parameters:

      if(nuclear.eq.'Fermi') then

       apar = a_fermi / fermi_in_cm * bohr_radius_in_cm
       cpar = c_fermi / fermi_in_cm * bohr_radius_in_cm

       if(a_fermi.eq.0.d0) then
        apar = t_fermi/(four*LOG(three))
        a_fermi = apar * fermi_in_cm / bohr_radius_in_cm
       end if

       if(c_fermi.eq.0.d0) then
        CALL GETCPR(rrms,apar,cpar)
        c_fermi = cpar * fermi_in_cm / bohr_radius_in_cm
       end if
      end if

      End Subroutine Read_nuclear


!======================================================================
      Real(8) Function Z_nuclear(r)
!======================================================================
!     Nucleus charge distribution
!----------------------------------------------------------------------
      Use DBS_nuclear, b => r_uniform, a => a_fermi, c => c_fermi, &
                       Z => atomic_number
      Implicit none
      Real(8), intent(in) :: r
      Real(8) :: ro

      if(r.le.zero) then
       Z_nuclear = zero
      elseif(nuclear.eq.'point') then       ! point nucleus
       Z_nuclear = zero                                          ! ???
       if(r.le.0.0219)  Z_nuclear = Z       ! after GRASP
      elseif(nuclear.eq.'uniform') then     ! uniform distribution
       ro_uniform = Z/(4*pi/3*b**3)
       Z_nuclear = zero
       if(r.le.b)  Z_nuclear = ro           ! ???
      elseif(nuclear.eq.'Fermi') then       ! Fermi distribution
       if(ro_fermi.eq.zero) Call  DEF_ro_fermi
       Z_nuclear = ro_fermi/(one+exp((r-c)/a))
      else
        Stop 'Z_nuclear: unknown charge distribution '
      end if

      End Function Z_nuclear


!======================================================================
      Real(8) Function V_nuclear(r)
!======================================================================
!     Nucleus potential at point "r" for different charge distributions
!     F. A. Parpia and A. K. Mohanty
!     "Relativistic basis set calculations for atoms with Fermi nuclei"
!     Phys Rev A46,3735 (1992)
!----------------------------------------------------------------------
      Use DBS_nuclear, Z => atomic_number, b => r_uniform, &
                       a => a_fermi, c => c_fermi

      Implicit none
      Real(8), intent(in) :: r
      Real(8) :: V, N, ra,rc, ac,ac2,ac3, s,s2,s3, p
      Real(8), external :: SKFUN

      if(r.le.zero) then

       V = zero

      elseif(nuclear.eq.'point') then       ! point nucleus

       V = -Z / r

      elseif(nuclear.eq.'uniform') then     ! uniform distribution

       if(r.gt.b) then
        V = -Z / r
       else
        V = -(three*Z)/(two*b)*(one - r*r/(three*b*b))
       end if

      elseif(nuclear.eq.'Fermi') then       ! Fermi distribution

        ac = a/c; ac2=ac*ac; ac3=ac2*ac
        p = pi*pi*ac2; s = six*ac3*SKFUN(3,-c/a)
        N = one + p - s
        rc = r/c

       if(r.le.c) then
        V = -Z / r
        ra = (r-c)/a; s2 = ac2*SKFUN(2,ra); s3 = ac3*SKFUN(3,ra)
        V = (three - rc*rc + p)/two - three*s2 - s/rc + six/rc*s3
        V = -(Z*V)/(N*c)
       else
        ra = (c-r)/a; s2 = ac2*SKFUN(2,ra); s3 = ac3*SKFUN(3,ra)
        V = N + six*s3 + three*rc*s2
        V = -Z/r * (V/N)
       end if

      else

        Stop 'V_nuclear: unknown charge distribution '

      end if

      V_nuclear = V

      End Function V_nuclear


!=======================================================================
      Real(8) Function SKFUN (k,x)
!=======================================================================
!              infinity   (-1)^n e^nx
!      S (x) =   Sum      -----------
!       k        n=1          n^k
!======================================================================
      Implicit none
      Integer, intent(in) :: k
      Real(8), intent(in) :: x
      Integer :: n
      Real(8) :: eps, base, dnum, dk, delta

      SKFUN = 0.d0
      if(x.lt.-500d0) Return

      eps = 10.d0*epsilon(1.d0)
      BASE = -exp(x)
      DNUM = BASE
      SKFUN = BASE
      n = 1
      Do
       n = n + 1
       DNUM = DNUM*BASE
       dk = n; dk = dk ** k
       DELTA = DNUM/dk
       SKFUN = SKFUN+DELTA
       if(abs(DELTA/SKFUN).le.eps) Exit
      End do

      End Function SKFUN


!=====================================================================
      Real(8) Function ESTRMS (APARM,CPARM)
!=====================================================================
!     Determines the root mean square radius for a Fermi nucleus with
!     given  the parameters `c' (CPARM) and `a' (APARM).
!     Based on F. A. Parpia and A. K. Mohanty  Phys Rev A 46,3735(1992)
!     and W. R. Johnson, Lectures on Atomic Physics (2007)
!     Call(s) to: SKFUN.
!---------------------------------------------------------------------
      Use zconst

      Implicit none
      Real(8), intent(in) :: APARM,CPARM
      Real(8) :: A,P,X,DN,DD
      Real(8), external :: SKFUN

      A =  APARM/CPARM
      P =  PI*A
      X = -CPARM/APARM
      DN = 1.d0 + 10.d0/3.d0*P**2 + 7.d0/3.d0* P**4  &
               - 120.d0*A**5 * SKFUN(5,X)
      DD = 1.d0 + P**2 - 6.d0*A**3 * SKFUN(3,X)
      ESTRMS = CPARM * SQRT(3.d0*DN/5.d0/DD)

      End Function ESTRMS


!=====================================================================
      SUBROUTINE GETCPR (RRMS,APARM,CPARM)
!=====================================================================
!     Determines the parameter `c' (CPARM) for a Fermi nucleus,  given
!     the root mean square radius (RRMS) and the parameter `a' (APARM).
!     Call(s) to: ESTRMS.
!----------------------------------------------------------------------
      Use zconst

      Implicit none
      Real(8), intent(in) :: RRMS,APARM
      Real(8), intent(out) :: CPARM
      Real(8) :: ACCY,CPMIN,CPMAX,CPTRY,RMSTRY
      Real(8), external :: ESTRMS

! ... accuracy parameter

      ACCY = two*epsilon(one)

! ... bracket CPARM with a lower and upper limits:

      CPMIN = RRMS
    1 CPMIN = half*CPMIN
      if (ESTRMS(APARM,CPMIN) .GT. RRMS) GOTO 1

      CPMAX = RRMS
    2 CPMAX = two*CPMAX
      if (ESTRMS(APARM,CPMAX) .LT. RRMS) GOTO 2

! ... find CPARM by the method of bisection

    3 CPTRY = half*(CPMAX+CPMIN)

      RMSTRY = ESTRMS(APARM,CPTRY)

      if (RMSTRY .GT. RRMS) then
        CPMAX = CPTRY
      else
        CPMIN = CPTRY
      END if

      if ( ( (CPMAX-CPMIN)/(CPMAX+CPMIN) .GT. ACCY )   .AND. &
           ( ABS(RMSTRY-RRMS)/RRMS       .GT. ACCY ) )  GOTO 3

      CPARM = CPTRY

      End Subroutine GETCPR

