!=======================================================================
      Subroutine CGAMMA (ARGR,ARGI,RESR,RESI)
!=======================================================================
!   This subroutine returns in (RESr,RESi) the complex Gamma function  !
!   of the complex argument (ARGr,ARGi).                               !
!                                                                      !
!   Only RESR is nonzero if ARGI is zero.                              !
!                                                                      !
!   The  ARCTAN function required must return angles (in radians) in   !
!   the range  [0,2*pi).                                               !
!                                                                      !
!   Written by Farid A Parpia, at Oxford    Last update: 06 Oct 1992   !
!=======================================================================

      Implicit real(8) (A-H,O-Z)

!----------------------------------------------------------------------
!   These are the Bernoulli numbers B02, B04, ..., B14, expressed as   
!   rational numbers. From Abramowitz and Stegun, p. 810.              

      Real(8), Save :: PI, HLNTPI, eps, MAXEXP, MINEXP
      Real(8), Dimension(7), SAVE :: FN,FD
      Data FN/ 1.0d0, -1.0d0,  1.0d0, -1.0d0,  5.0d0, -691.0d0, 7.0d0/
      Data FD/ 6.0d0, 30.0d0, 42.0d0, 30.0d0, 66.0d0, 2730.0d0, 6.0d0/

!----------------------------------------------------------------------
      Logical, save :: FIRST, NEGARG
      Data FIRST/.TRUE./
      Integer :: i

!   On the first entry to this routine, set up the constants required
!   for the reflection formula (cf. Abramowitz and Stegun 6.1.17) and
!   Stirling's approximation (cf. Abramowitz and Stegun 6.1.40).

      IF(FIRST) THEN

         PI = ACOS(-1.d0)
         HLNTPI = 0.5D0*LOG (PI+PI)
         eps = 2*EPSILON(1.d0)
         MAXEXP = MAXEXPONENT(1.d0) * LOG(2.d0)
         MINEXP = MINEXPONENT(1.d0) * LOG(2.d0)

         Do i = 1,7;  DI = i+i
           FN(I) = FN(I)/(FD(I)*DI*(DI-1.0d0))
         End do

         FIRST = .FALSE.

      END IF

!----------------------------------------------------------------------
! ... Cases where the argument is real

      IF (ARGI .EQ. 0.d0) THEN

! ... Cases where the argument is real and negative

         IF (ARGR .LE. 0.d0) THEN

! ... Stop with an error message if the argument is too near a pole

            IF (ABS (DBLE (NINT (ARGR))-ARGR) .LE. eps) &
               STOP 'CGAMMA: Argument too close to a pole.'
            
!   Otherwise use the reflection formula (Abramowitz and Stegun 6.1.17)
!   to ensure that the argument is suitable for Stirling's formula

            ARGUM = PI/(-ARGR*SIN(PI*ARGR))
            CLNGI = 0.d0
            IF (ARGUM .LT. 0.d0) THEN
              ARGUM = -ARGUM
              CLNGI = PI
            ENDIF
            FACNEG = LOG (ARGUM)
            ARGUR = -ARGR
            NEGARG = .TRUE.

! ...  Cases where the argument is real and positive

         ELSE

            CLNGI = 0.0d0
            ARGUR = ARGR
            NEGARG = .FALSE.

         END IF

!   Use Abramowitz and Stegun formula 6.1.15 to ensure that
!   the argument in Stirling's formula is greater than 10

         OVLFAC = 1.0d0
         Do
          IF (ARGUR .GE. 10.0d0) Exit
          OVLFAC = OVLFAC*ARGUR
          ARGUR = ARGUR+1.0d0
         End do

!   Now use Stirling's formula to compute Log (Gamma (ARGUM))

         CLNGR = (ARGUR-0.5d0)*LOG(ARGUR)-ARGUR+HLNTPI
         FAC = ARGUR
         OBASQ = 1.0d0/(ARGUR*ARGUR)
         Do I = 1,7
            FAC = FAC*OBASQ; CLNGR = CLNGR+FN(I)*FAC
         End do

!   Include the contributions from the recurrence and reflection
!   formulae

         CLNGR = CLNGR-LOG (OVLFAC)
         IF (NEGARG) CLNGR = FACNEG-CLNGR

      ELSE

!   Cases where the argument is complex

         ARGUR = ARGR
         ARGUI = ARGI
         ARGUI2 = ARGUI*ARGUI

!   Use the recurrence formula (Abramowitz and Stegun 6.1.15)
!   to ensure that the magnitude of the argument in Stirling's
!   formula is greater than 10

        OVLFR = 1.0d0
        OVLFI = 0.0d0
        Do
         ARGUM = SQRT (ARGUR*ARGUR+ARGUI2)
         IF (ARGUM .GE. 10.d0) Exit
         TERMR = OVLFR*ARGUR-OVLFI*ARGUI
         TERMI = OVLFR*ARGUI+OVLFI*ARGUR
         OVLFR = TERMR
         OVLFI = TERMI
         ARGUR = ARGUR+1.0d0
        End do

!   Now use Stirling's formula to compute Log (Gamma (ARGUM))

         ARGUR2 = ARGUR*ARGUR
         TERMR = 0.5d0*LOG (ARGUR2+ARGUI2)
         TERMI = ATAN2 (ARGUI,ARGUR)
         if(TERMI.lt.0.d0) TERMI = TERMI + PI + PI 
         CLNGR = (ARGUR-0.5d0)*TERMR - ARGUI*TERMI-ARGUR+HLNTPI
         CLNGI = (ARGUR-0.5d0)*TERMI + ARGUI*TERMR-ARGUI
         FAC = (ARGUR2+ARGUI2)**(-2)
         OBASQR = (ARGUR2-ARGUI2)*FAC
         OBASQI = -2.0d0*ARGUR*ARGUI*FAC
         ZFACR = ARGUR
         ZFACI = ARGUI
         Do I = 1,7
            TERMR = ZFACR*OBASQR-ZFACI*OBASQI
            TERMI = ZFACR*OBASQI+ZFACI*OBASQR
            FAC = FN(I)
            CLNGR = CLNGR+TERMR*FAC
            CLNGI = CLNGI+TERMI*FAC
            ZFACR = TERMR
            ZFACI = TERMI
         End do

!   Add in the relevant pieces from the recurrence formula

         CLNGR = CLNGR - 0.5d0 * LOG (OVLFR*OVLFR+OVLFI*OVLFI)
         A = ATAN2(OVLFI,OVLFR); if(A.lt.0.d0) A = A + PI + PI
         CLNGI = CLNGI - ATAN2(OVLFI,OVLFR)
 
      END IF
 
!   Now exponentiate the complex Log Gamma function to get
!   the complex Gamma function
 
      IF (CLNGR.ge.MAXEXP.or.CLNGR.le.MINEXP) then
          write(*,*) CLNGR,MAXEXP,MINEXP
         STOP 'CGAMMA: Argument to exponential function out of range.'
      End if

      FAC = EXP(CLNGR)
      RESR = FAC*COS(CLNGI)
      RESI = FAC*SIN(CLNGI)

      End Subroutine CGAMMA


!======================================================================
      Function RGAMMA(X)
!======================================================================
!     COMPUTES GAMMA(X) TO ABOUT 14 FIGURES
!     presented by Peter J. Mohr
!----------------------------------------------------------------------
      Use zconst

      Implicit none
      Real(dp), intent(in) :: X
      Real(dp) :: RGAMMA, Y,Z, A1,A2,A3,A4,A5,A6,A7, C,D,G

      Y=X 
      G=one
    1 IF(Y.GT.seven)  go to 2
      G=G*Y
      Y=Y+one
      go to 1
    2 Z=-one/(Y*Y)

      C =.398942280401433e+0_dp
      D =.295506535947712e-1_dp

      A1=.282001658833287e+1_dp
      A2=.940005529444291e-1_dp
      A3=.268573008412655e-1_dp
      A4=.201429756309491e-1_dp
      A5=.284850160437664e-1_dp
      A6=.648894925920085e-1_dp
      A7=.216924352948683e+0_dp

      RGAMMA = (Y-half)*LOG(Y) - Y - LOG(C*G) + &
               (A1+Z*(A2+Z*(A3+Z*(A4+Z*(A5+Z*(A6+Z*(A7+Z)))))))*D/Y

      RGAMMA = EXP(RGAMMA)

      End Function RGAMMA
