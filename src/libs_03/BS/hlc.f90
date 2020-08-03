!======================================================================
      DOUBLE PRECISION FUNCTION HLC(I,J)
!======================================================================
!     COMPUTES HL(I,J) MODIFIED BY THE INTERACTIONS WITH CLOSED SHELLS
!----------------------------------------------------------------------
      USE spline_orbitals, l => lbs
      USE spline_atomic, nclosd => kclosd

      IMPLICIT NONE
      INTEGER, INTENT(in) :: i,j
      INTEGER :: ip, k, kmin, kmax
      REAL(8) :: CA, CB
      REAL(8), EXTERNAL :: ZCB, RKy, BHL

      HLC = BHL(I,J)

! ... direct interaction:

      DO IP = 1,NCLOSD
        CA = -2.d0*(4*L(IP)+2)
        HLC = HLC + CA*RKY(I,IP,J,IP,0)
      END DO

! ... exchange contribution:

      kmin = 1000;  kmax = 0
      DO IP = 1,NCLOSD
       k = iabs(L(I)-L(IP));  if(k.lt.kmin) kmin = k
       k =      L(I)+L(IP) ;  if(k.gt.kmax) kmax = k
      End do

      DO k = kmin,kmax
       DO IP = 1,NCLOSD
         CB =  ZCB(L(I),K,L(IP)) * (4*L(IP)+2)
         if(CB.eq.0.d0) Cycle
         HLC = HLC + CB*RKY(I,IP,IP,J,k)
       END DO
      END DO

      END function HLC
