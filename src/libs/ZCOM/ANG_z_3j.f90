!====================================================================
      Real(8) Function Z_3j (j1,m1,j2,m2,j3,m3) 
!====================================================================
!     determines the value of the 3j-symbols without direct using of
!     factorials. The following expression for the 3j-symbols is used:
!         (A.P.JUCYS, A.A.BANDZAITIS, 1977)
!
!     3j{j1,m1,j2,m2,j3,m3} = delta(m1+m2,m3) * (2j3+1)^1/2 * {j1,j2,j3} *
!       sqrt[ (j1+m1)!*(j1-m1)!*(j2+m2)!*(j2-m2)!*(j3+m3)!*(j3-m3)! ]
!                         SUM(z) {   (-1)^z  /
!          [ z! *  (j1+j2-j3-z)! * (j1-m1-z)! * (j2-m2-z)! *
!                  (j3-j2+m1+z)! * (j3-j1-m2+z)! ] }
!
!     where {a,b,c}=sqrt[ (a+b-c)! * (a-b+c)! * (b-a+c)! / (a+b+c+1)! ]
!
!     If we introduce the auxiliary values a(i) and b(i)
!     (see below the text of program) then
!
!     3j =         (-1) ^ Sum[a(i)]
!          sqrt{ Pr[ (b(j)-a(i))! ] / [ Sum (b(j)-a(i))+1 ] }
!                  Sum(z) { (-1)^z  /
!          [  Pr[ (z-a(i))! ]  * Pr[ (b(j)-z)! ]   ] }
!
!     (below the moments are used in (2J+1)-representation)
!
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j1,m1,j2,m2,j3,m3
      Integer :: i,i_max,k,kk,m,iz,iz_min,iz_max
      Real(8) :: x,y,z
      Integer :: a(3),b(3),J(16) 

      Z_3j=0.0

      IF(M1+M2+M3-3.ne.0) Return    ! check of conservation rules
      J(1)= J1+J2-J3-1
      J(2)= J1-J2+J3-1
      J(3)= J2-J1+J3-1
      J(4)= J1+M1-2
      J(5)= J1-M1
      J(6)= J2-M2
      J(7)= J2+M2-2
      J(8)= J3+M3-2
      J(9)= J3-M3
      Do I=1,9
       IF(J(i).lt.0.or.mod(J(i),2).eq.1) Return
      End do

      a(1) = 0                         ! auxiliary values
      a(2) = (j2-j3-m1+1)/2
      a(3) = (j1-j3+m2-1)/2
      b(1) = (j1+j2-j3-1)/2
      b(2) = (j1-m1)/2
      b(3) = (j2+m2-2)/2

      IZ_min=MAX0(a(1),a(2),a(3))      ! limits of the sum
      IZ_max=MIN0(b(1),b(2),b(3))
      IF(IZ_max.LT.IZ_min) Return

      Do I=1,3                         ! constant factorial parameters
      Do K=1,3
       J(I+3*K-3)=b(i)-a(k)
      End do
      End do
      J(10)=(j1+j2+j3-3)/2+1

      Do I=1,3
       J(I+10)=IZ_min-a(i)               ! initial factorial parameters
       J(I+13)=b(i)-IZ_min               ! in the sum
      End do

      Z=0.0
      DO IZ=IZ_min,IZ_max                 ! summation

       I_max=0                            ! max. factorial
       Do I=1,16
        if(J(i).gt.I_max) I_max=J(i)
       End do

       Y=1.0
       DO I=2,I_max         ! estimation of one term in sum
        K=0                 ! K - the extent of the integer I in term
        DO M=1,9
         IF(J(M).GE.I) K=K+1
        End do
        IF(J(10).GE.I) K=K-1
        DO M=11,16
         IF(J(M).GE.I) K=K-2
        End do
        IF(K.EQ.0) Cycle

        X=DBLE(I)                   ! Y = Y * I ** K/2
        KK=IABS(K)/2
        IF(KK.GT.0) THEN
         DO M=1,KK
          IF(K.GT.0) Y=Y*X
          IF(K.LT.0) Y=Y/X
         END DO
        END IF
        IF(mod(K,2).EQ.+1) Y=Y*SQRT(X)
        IF(mod(K,2).EQ.-1) Y=Y/SQRT(X)
       End do

       IF(mod(IZ,2).eq.1) Y=-Y
       Z=Z+Y

       Do I=11,13                  ! new factorial parameters in sum
        J(I)=J(I)+1
       End do
       DO I=14,16
        J(I)=J(I)-1
       End do

      End do                       ! end of summation

      K=a(1)+a(2)+a(3)
      if(mod(k,2).ne.0) Z=-Z
      Z_3j=Z

      END Function Z_3j


!====================================================================
      Real(8) Function Z_3jj(j1,m1,j2,m2,j3,m3)
!====================================================================
!     3j-symbols for L-moments (without 2J+1 representation
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j1,m1,j2,m2,j3,m3
      Real(8), external :: Z_3j
       
      Z_3jj=Z_3j(j1+j1+1,m1+m1+1,j2+j2+1,m2+m2+1,j3+j3+1,m3+m3+1)

      End Function Z_3jj


!====================================================================
      Real(8) Function Z_3j2(j1,m1,j2,m2,j3,m3)
!====================================================================
!     3j-symbols for moments in 2J representation
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j1,m1,j2,m2,j3,m3
      Real(8), External :: Z_3j
       
      Z_3j2=Z_3j(j1+1,m1+1,j2+1,m2+1,j3+1,m3+1)

      End Function Z_3j2


!====================================================================
      Real(8) Function CLEBSH(J1,M1,J2,M2,J,M)
!====================================================================
!     determines the Clebsh-Gordon coefficients through the 3j-symbol
!     (the moments are used in (2J+1)-representation)
!
!     Call:  Z_3j
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: J, M, J1, M1, J2, M2
      Real(8), external :: Z_3j

      Clebsh=(-1)**((j1-j2+m-1)/2)*sqrt(DBLE(J))*   &
             Z_3j(j1,m1,j2,m2,J,-m+2)

      END Function CLEBSH


!======================================================================
      Real(8) Function CLEBCH(L1,M1,L2,M2,L,M)
!======================================================================
!     determines the Clebsh-Gordon coefficients through the 3j-symbol
!     (the moments are in L-representation)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: L, M, L1, M1, L2, M2 
      Real(8), External :: Clebsh

      Clebch = Clebsh(l1+l1+1,m1+m1+1,l2+l2+1,m2+m2+1,l+l+1,m+m+1)

      End Function CLEBCH


!======================================================================
      Real(8) Function CLEBSH2(L1,M1,L2,M2,L,M)
!======================================================================
!     determines the Clebsh-Gordon coefficients through the 3j-symbol
!     (the moments are in 2J-representation)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: L, M, L1, M1, L2, M2 
      Real(8), External :: Clebsh

      Clebsh2 = Clebsh(l1+1,m1+1,l2+1,m2+1,l+1,m+1)

      End Function CLEBSH2


!======================================================================
      Real(8) Function Z_3j0 (K1,K2,K3)
!======================================================================
!     { k1 k2 k3}   -  3j-symbol with zero M-values
!     {  0  0  0}
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Integer, intent(in) :: k1,k2,k3

      Z_3j0=0.d0

      M=K1+K2+K3; N=M/2; IF(N+N.NE.M) RETURN; M=M+1

      N1=N-K1; N2=N-K2; N3=N-K3
      IF(N1.LT.0.OR.N2.LT.0.OR.N3.LT.0) RETURN
      M1=N1+N1; M2=N2+N2; M3=N3+N3

      A=1.d0
      DO I=2,M

      K=-1                
      IF(M1.GE.I) K=K+1; IF(M2.GE.I) K=K+1; IF(M3.GE.I) K=K+1
      IF(N .GE.I) K=K+2
      IF(N1.GE.I) K=K-2; IF(N2.GE.I) K=K-2; IF(N3.GE.I) K=K-2

      IK=IABS(K); B=DBLE(I)
      IF(K.GT.0) THEN
       DO J=1,IK; A=A*B; END DO
      ELSEIF(K.LT.0) THEN
       DO J=1,IK; A=A/B; END DO
      END IF

      END DO

      Z_3j0=SQRT(A); IF(mod(N,2).eq.1) Z_3j0=-Z_3j0

      End Function Z_3j0

!======================================================================
      Real(8) Function ZCB(K1,K2,K3)
!======================================================================
!     CB =  3j(k1,0,k2,0,k3,0)**2
!
!     3j =  sqrt[ (2g-2k1)! (2g-2k2)! (2g-2k3)! / (2g+1)! ] *
!                   g! / [ (g-k1)! (g-k2)! (g-k3)! ],
!
!     where 2g = k1+k2+k3
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k1,k2,k3
      Real(8) :: A,B,C
      Integer :: I,J, K,L, M,N, N1,N2,N3, M1,M2,M3

      ZCB=0.0;  M=K1+K2+K3;  N=M/2; if(N+N.NE.M) Return
      N1=N-K1; N2=N-K2; N3=N-K3
      IF(N1.LT.0.OR.N2.LT.0.OR.N3.LT.0) Return

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


!======================================================================
      Real(8) Function CLB(K1,K2,K3)
!======================================================================
!     Clebsh(k1,0,k2,0;k3,0) = (-1)^(k1-k2) [k3]^1/2 3j(k1,0,k2,0;k3,0) 
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k1,k2,k3
      Real(8), external :: ZCB
      Real(8) :: A,B,C
      Integer :: I,J, K,L, M,N, N1,N2,N3, M1,M2,M3

      CLB = ZCB(k1,k2,k3)* (k3+k3+1)
      CLB = sqrt(CLB) * (-1)**(k1-k2)       ! sign?

      END FUNCTION CLB
