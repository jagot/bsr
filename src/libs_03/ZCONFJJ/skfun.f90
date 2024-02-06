! some programs from GRASP, used for nuclear parameters

!=====================================================================
      Real(8) Function ESTRMS (APARM,CPARM)
!=====================================================================
!  Determines the root mean square radius for a Fermi nucleus given
!  the parameters `c' (CPARM) and `a' (APARM).
!  Based on F. A. Parpia and A. K. Mohanty  Phys Rev A 46,3735(1992)
!  and W. R. Johnson, Lectures on Atomic Physics (2007)
!  Call(s) to: SKFUN.
!---------------------------------------------------------------------
      Use zconst, only: PI

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
      Subroutine GETCPR (RRMS,APARM,CPARM)
!=====================================================================
!  Determines the parameter `c' (CPARM) for a Fermi nucleus,  given
!  the root mean square radius (RRMS) and the parameter `a' (APARM).
!  Call(s) to: ESTRMS.
!----------------------------------------------------------------------
      Implicit none
      Real(8), intent(in) :: RRMS,APARM
      Real(8), intent(out) :: CPARM
      Real(8) :: ACCY,CPMIN,CPMAX,CPTRY,RMSTRY
      Real(8) :: one = 1.d0, two = 2.d0, half = 0.5d0
      Real(8), external :: ESTRMS

!     Accuracy parameter

      ACCY = two*epsilon(one)

!     Bracket CPARM with a lower and upper limits:

      CPMIN = RRMS
    1 CPMIN = half*CPMIN
      IF (ESTRMS(APARM,CPMIN) .GT. RRMS) GOTO 1

      CPMAX = RRMS
    2 CPMAX = two*CPMAX
      IF (ESTRMS(APARM,CPMAX) .LT. RRMS) GOTO 2

! ... find CPARM by the method of bisection

    3 CPTRY = half*(CPMAX+CPMIN)

      RMSTRY = ESTRMS(APARM,CPTRY)

      IF (RMSTRY .GT. RRMS) THEN
         CPMAX = CPTRY
      ELSE
         CPMIN = CPTRY
      END IF

      IF ( ( (CPMAX-CPMIN)/(CPMAX+CPMIN) .GT. ACCY )   .AND. &
           ( ABS(RMSTRY-RRMS)/RRMS       .GT. ACCY ) )  GOTO 3

      CPARM = CPTRY

      End Subroutine GETCPR


!=======================================================================
      Real(8) Function SKFUN (k,x)
!=======================================================================
!                                    n  nx
!                     infinity   (-1)  e
!             S (x) =   Sum      -------
!              k        n=1          k
!                                   n
!======================================================================
      Implicit none
      Integer, intent(in) :: k
      Real(8), intent(in) :: x
      Integer :: n
      Real(8) :: eps, base, dnum, dk, delta

      eps = 10.d0*epsilon(1.d0)
      BASE = -exp(x)
      DNUM = BASE
      SKFUN = BASE
      if(BASE.eq.0.d0) Return

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
