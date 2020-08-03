!======================================================================
      Real(8) Function ZCLKL(K1,K2,K3)
!======================================================================
!     reduced matrix elements of spherical garmonics:
!
!     <k1||C[k2]||k3> = (-1)^k1 * sqtr[ (2*k1+1) * (2*k3+1) ] *
!                     3j(k1,0,k2,0,k3,0)
!
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(In) :: k1,k2,k3
      Integer :: M, N
      Real(8), external :: ZCB

      M=K1+K2+K3
      N=M/2
      IF(N+N.NE.M) Return

      ZCLKL = ZCB(K1,K2,K3)
      ZCLKL = SQRT(ZCLKL*(K1+K1+1)*(K3+K3+1))
      IF(mod(N+K1,2).eq.1) ZCLKL=-ZCLKL

      End Function ZCLKL
