!=====================================================================
      REAL(8) FUNCTION ZINT(icase,k,i1,i2,i3,i4)
!=====================================================================
!
!     Calculation of radial integral according its type  'icase'
!
!     This function serves only for connection with Breit_no, or
!     BSR_breit, output in INT_BNK.
!---------------------------------------------------------------------

      USE RADIAL

      IMPLICIT REAL(8) (A-H,O-Z)

      Real(8), External :: RK,NK,MK,TK,WK,HLC,ZETA

      Integer(4), Intent(in) :: icase,k,i1,i2,i3,i4

      if(icase.le.0.or.icase.gt.10) Stop ' ZINT: icase is out of range'

      Select case (icase)

      Case(1);  ZINT = RK(i1,i2,i1,i2,k)       !  FK

      Case(2);  ZINT = RK(i1,i2,i2,i1,k)       !  GK

      Case(3);  ZINT = TK(I1,I2,I3,I4,K)       !  TK

      Case(4);  ZINT = MK(I1,I2,I3,I4,K)       !  MK

      Case(5);  ZINT = RK(I1,I2,I3,I4,K)       !  RK

      Case(6);  ZINT = HLC(I1,I3)              !  HL

      Case(7);  ZINT = ZETA(I1,I3)             !  ZETA

      Case(8);  ZINT = NK(I1,I2,I3,I4,K)       !  NK

      Case(9);  ZINT = VK(I1,I2,I3,I4,K)       !  VK                     !  VK

      Case(10); ZINT = NK(I1,I2,I3,I4,K)       !  NK

      End select

      End Function Zint



!----------------------------------------------------------------------
      REAL(8) FUNCTION MK (I,J,II,JJ,K)
!----------------------------------------------------------------------
!
!     returns value of MK-integral:
!
!     Mk(I,J,II,JJ) = (Nk(I,J,II,JJ) + Nk(J,I,JJ,II))/2
!
!----------------------------------------------------------------------

      USE RADIAL

      Implicit none

      Integer(4), Intent(in) :: i,ii,j,jj,k

      Real(8), External :: NK

      MK = (NK(I,J,II,JJ,K) + NK(J,I,JJ,II,K))/2.d0

      END FUNCTION MK



!=====================================================================
      REAL(8) FUNCTION WK (I,J,II,JJ,K)
!=====================================================================
!
!     Evaluate WK integral:
!
!     W(I,J,II,JJ,K) =     2 * V(I,J,II,JJ,K)
!                     -(k+1) * N(I,J,II,JJ,K+1)
!                     +(k+2) * N(J,I,JJ,II,K-1)
!--------------------------------------------------------------------

      IMPLICIT NONE

      Integer(4), Intent(in) :: I,II,J,JJ,K
      Real(8), External :: VK,NK

      WK = 2 * VK(I,J,II,JJ,K) -(k+1) * NK(I,J,II,JJ,K+1) &
                               +(k+2) * NK(J,I,JJ,II,K-1)

      END FUNCTION WK
