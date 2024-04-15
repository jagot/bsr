!======================================================================
      SUBROUTINE INITR
!======================================================================
!
!     set the exponential radial grid and related variables
!
!----------------------------------------------------------------------

      USE RADIAL

      IMPLICIT NONE
      Integer(4) :: i

      if(.not.allocated(R)) Allocate(R(NR),RR(NR),R2(NR),ZK(NR),YK(NR))

      if(Z.le.0.d0) Stop 'INITR: Z <= 0 '

      H1 = H*D2/D3
      EH = DEXP(-H)

      DO I=1,NR
       R(I)=EXP(RHO+(I-1)*H)/Z
       RR(I)=R(I)*R(I)
       R2(I)=SQRT(R(I))
      END DO


      END SUBROUTINE INITR