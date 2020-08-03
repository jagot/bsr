!---------------------------------------------------------------------
      Real(8) FUNCTION GRAD(I,J)
!---------------------------------------------------------------------
!
!     returns the dipole integral between i-th and j-th orbitals
!     in the V-form:  INT{ P_i(r) [d/dr +- l>/r] P_j(r); dr}
!                     li > lj,  insert (-), then li < lj
!---------------------------------------------------------------------

      Use RADIAL, L => lro

      Implicit real(8) (A-H,O-Z)
      Integer(4), Intent(in) :: I,J

      if(IABS(L(I)-L(J)).ne.1) STOP  ' Grad: LI - LJ <> 1'
      LL = MAX0(L(I),L(J))

      II = I; JJ = J
      if (L(I).lt.L(J)) then;  II = J; JJ = I; end if

      K = 1
      A1 = (LL + D5)/(LL*(LL+1)*(2*LL+1))
      GRAD = R(1)*P(1,I)*P(1,J) * (D1+A1*Z*R(1))
      DL = D5*P(1,I)*P(1,J)*R(1)

      K = 2
      F1 = D5*(P(K+1,II) - P(K-1,II))
      F2 = P(K+1,II) - D2*P(K,II) + P(K-1,II)

      G0 = P(K,JJ)*R(K)
      G1 = D5*(P(K+1,JJ)*R(K+1) - P(K-1,JJ)*R(K-1))
      G2 = P(K+1,JJ)*R(K+1) - D2*P(K,JJ)*R(K) + P(K-1,JJ)*R(K-1)

      GRAD = GRAD + D2*F1*G0 +(D2*F2*G1 + F1*G2)/D3

      DL = DL + D2*P(K,II)*P(K,JJ)*R(K) + P(K+1,II)*P(K+1,JJ)*R(K+1)

      MM = MIN0(mro(I)+1,mro(J)+1,NR-2)

      DO K = 4,MM,2

       F1 = D5*(P(K+1,II) - P(K-1,II))
       F2 = P(K+1,II) - D2*P(K,II) + P(K-1,II)
       F3 = D5*(P(K+2,II) - P(K-2,II)) - D2*F1
       F4 = P(K+2,II) + P(K-2,II) - D4*(P(K+1,II) + P(K-1,II)) &
          + D6*P(K,II)

       G0 = P(K,JJ)*R(K)
       G1 = D5*(P(K+1,JJ)*R(K+1) - P(K-1,JJ)*R(K-1))
       G2 = P(K+1,JJ)*R(K+1) - D2*P(K,JJ)*R(K) + P(K-1,JJ)*R(K-1)
       G3 = D5*(P(K+2,JJ)*R(K+2) - P(K-2,JJ)*R(K-2)) - D2*G1
       G4 = P(K+2,JJ)*R(K+2) + P(K-2,JJ)*R(K-2) - D4*(P(K+1,JJ)*R(K+1) &
          + P(K-1,JJ)*R(K-1)) + D6*P(K,JJ)*R(K)

       GRAD = GRAD + D2*F1*G0 +(D2*F2*G1 + F1*G2)/D3 &
            - (F1*G4-F4*G1 + D4*(F2*G3-F3*G2))/90.D0

       DL = DL + D2*P(K,II)*P(K,JJ)*R(K) + P(K+1,II)*P(K+1,JJ)*R(K+1)

      End do

      GRAD = GRAD + (LL+D5)*DL*H1

      IF (II .EQ. J) GRAD = - GRAD  

      END FUNCTION GRAD


!---------------------------------------------------------------------
      Real(8) FUNCTION ZGRAD(I,J)
!---------------------------------------------------------------------
!
!     returns the dipole integral between i-th and j-th orbitals
!     in the V-form:  INT{ P_j(r) [d/dr + l>/r] P_i(r); dr}
!                     li > lj,  insert (-), then li < lj
!                  
!---------------------------------------------------------------------

      Use RADIAL, L => lro

      Implicit real(8) (A-H,O-Z)
      Integer(4), Intent(in) :: I,J

      if(IABS(L(I)-L(J)).ne.1) STOP  ' ZGRAD: LI - LJ <> 1'
      LL = MAX0(L(I),L(J));  if(L(i).gt.L(j)) LL=-LL

      IP = nrf+1;  Call P_derive(J,IP)
	  
      G = QUADR(I,IP,0); U = QUADR(I,J,-1)

      ZGRAD = G + (LL+1)*U
 
      END FUNCTION ZGRAD

