!======================================================================
      Real(8) Function E_abc(a,b,c)
!======================================================================
!
!     E(a,b,c) = (a+b+c+1)!((a+b-c)/2)!((a-b+c)/2)!((-a+b+c)/2)!
!                /(a+b-c)!   /(a+b+c)!
!
!     Parameters a,b,c are restricted by triangle relation {a,b,c} 
!     and a+b+c must be even
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: a,b,c
      Real(8) :: A,B,C
      Integer :: I,J, K,L, M, N, N1,N2,N3, M1,M2

      E_abc=0.0; if(mod(a+b+c,2).ne.0) Return 
      if(itra(a,b,c).eq.0) Return
      if(a.lt.0.or.b.lt.0.or.c.lt.0) Return

      N  =  a + b + c;       M  = N/2  
      N1 = (a + b - c)/2;    M1 = N1*2 
      N2 = (a - b + c)/2 
      N3 = (b - a + c)/2 


      A = N + 1
      DO I = 1,M
       K = +1
       IF(N1.GE.I) K=K+1
       IF(N2.GE.I) K=K+1
       IF(N3.GE.I) K=K+1
       IF(M .GE.I) K=K-1
       IF(M1.GE.I) K=K-1

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

