!======================================================================
      Integer Function Iopen(n,E,ETARG)
!======================================================================
!     number of open target states for given energy E
!----------------------------------------------------------------------
      Implicit none
      Integer,intent(in) :: n
      Real(8),intent(in) :: E
      Real(8),intent(in) :: Etarg(n)
      Integer :: i

      Iopen=0;  Do i=1,n; if(E.gt.ETARG(i)) Iopen=i; End do

      End Function Iopen


!=======================================================================
      Subroutine ZOPEN(ETOT,nast,ENAT,NCONAT,nch,nopen,EN)
!=======================================================================
!      DETERMINE NUMBER OF OPEN CHANNELS AND CHANNELS ENERGIES
!      INPUT:
!       ETOT       = ELECTRON ENERGY IN RYDS
!       nast       = NUMBER OF TARGET STATES
!       ENAT       = TARGET ENERGIES IN RYDS
!       NCONAT     = NUMBER OF CHANNELS COUPLED TO EACH TARGET STATE
!       NCH        = NUMBER OF CHANNELS
!      RETURNS:
!       NOPEN      = NUMBER OF OPEN CHANNELS
!       EN         = CHANNEL ENERGIES IN RYDS.
!-----------------------------------------------------------------------
      Implicit none
      Integer, intent(in)  :: nast, nch
      Integer, intent(out) :: nopen
      Integer, intent(in)  :: NCONAT(nast)
      Real(8), intent(in)  :: ETOT
      Real(8), intent(in)  :: ENAT(nast)
      Real(8), intent(out) :: EN(nch)
      Integer :: K,N,NC
      Real(8) :: ECH

      NOPEN = 0
      K = 0
      DO N = 1,nast
        IF (NCONAT(N).EQ.0) Cycle
        ECH = ETOT - ENAT(N)
        IF (ECH.GT.0.d0) NOPEN = NOPEN + NCONAT(N)
        DO NC = 1,NCONAT(N)
          K = K + 1
          EN(K) = ECH
        END DO
      END DO

      End Subroutine ZOPEN

