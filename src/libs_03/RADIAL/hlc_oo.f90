!======================================================================
      REAL(8) FUNCTION HLC_oo (I,J)
!======================================================================
!
!     COMPUTES MODIFICATION OF HL(I,J) DUE TO ORBIT-ORBIT INTERACTION
!     WITH THE CLOSED SHELLS
!
!----------------------------------------------------------------------

      USE RADIAL, l => lro

      IMPLICIT NONE
      INTEGER(4),INTENT(in) :: i,j
      Integer(4) :: k,ip,LI,LP
      Real(8) :: C,CB,CN
      Real(8), External :: RK_oo, NK, ZCB

      HLC_oo = D0

      LI=L(i)
      DO IP = 1,KCLOSD
       LP=L(IP)
       C = D0
       DO K = IABS(LI-LP),LI+LP,2
        if(k.eq.0) Cycle
        CB = -ZCB(LI,K,LP)/2
        if(CB.eq.D0) Cycle
        C = C + CB*RK_oo(I,IP,IP,J,K)
        if(LI.eq.0.or.LP.eq.0) Cycle
        CN = CB * (LI+LP+k+2)*(LI+LP-k)*(LI-LP+k+1)*(LP-LI+k+1) &
                / ((k+1)*(k+2))
        C = C + CN*(NK(I,IP,IP,J,K) + NK(IP,I,J,IP,K))/D2   
       END DO
       HLC_oo = HLC_oo - 2*(4*LP+2)*C
      END DO     

      END FUNCTION HLC_oo
