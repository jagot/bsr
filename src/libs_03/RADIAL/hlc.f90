!======================================================================
      REAL(8) FUNCTION HLC (I,J)
!======================================================================
!
!     COMPUTES MODIFICATION OF HL(I,J) DUE TO INTERACTION
!     WITH THE CLOSED SHELLS
!
!----------------------------------------------------------------------

      USE RADIAL, L => lro

      IMPLICIT NONE
      INTEGER(4),INTENT(in) :: i,j
      Integer(4) :: k,LI,LP,ip
      Real(8) :: C
      Real(8), External :: RK, ZCB, HL, HLC_oo

      HLC = HL(i,j)

      LI=L(i)
      DO IP = 1,KCLOSD
       LP=L(IP)

!  ... direct interaction:

       C = RK(I,IP,J,IP,0)

!  ... exchange interaction:

       DO K = IABS(LI-LP),LI+LP,2
        C = C -ZCB(LI,K,LP)*RK(I,IP,IP,J,K)/D2
       End do

       HLC = HLC - 2*(4*LP+2)*C

      END DO

! ... correction for orbit-orbit interaction:

      if(oo) HLC = HLC + HLC_oo(I,J)

      END FUNCTION HLC
