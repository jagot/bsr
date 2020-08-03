!======================================================================
      Real(8) FUNCTION ZCB(K1,K2,K3)
!======================================================================
!
!     CB =  3j(k1,0,k2,0,k3,0)**2
!
!     3j =  sqrt[ (2g-2k1)! (2g-2k2)! (2g-2k3)! / (2g+1)! ] *
!                   g! / [ (g-k1)! (g-k2)! (g-k3)! ],
!
!     where 2g = k1+k2+k3
!----------------------------------------------------------------------

      Implicit none
      
      Integer(4), intent(in) :: k1,k2,k3
 
      Real(8) :: A,B,C
      Integer(4) :: I,J, K,L, M,N, N1,N2,N3, M1,M2,M3

      ZCB=0.0;  M=K1+K2+K3;  N=M/2; IF(N+N.NE.M) RETURN

      N1=N-K1; N2=N-K2; N3=N-K3

      IF(N1.LT.0.OR.N2.LT.0.OR.N3.LT.0) RETURN

      M1=N1+N1; M2=N2+N2; M3=N3+N3

      A = 1.d0/(M + 1)
      DO I = 1,N
       K = +1
       IF(M1.GE.I) K=K+1
       IF(M2.GE.I) K=K+1
       IF(M3.GE.I) K=K+1
       IF(N1.GE.I) K=K-2
       IF(N2.GE.I) K=K-2
       IF(N3.GE.I) K=K-2

       J = I + N
       L = -1 
       IF(M1.GE.J) L=L+1 
       IF(M2.GE.J) L=L+1 
       IF(M3.GE.J) L=L+1 

       B=I;  C=J;  A = A * B**K * C**L

      End do

      ZCB = A

      END FUNCTION ZCB


