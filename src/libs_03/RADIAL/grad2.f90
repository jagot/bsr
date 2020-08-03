!---------------------------------------------------------------------
      Real(8) FUNCTION GRAD2(I,J)
!---------------------------------------------------------------------
!
!     computes the dipole integral between i-th and j-th orbitals
!     in the V-form:  INT{ P_j(r) [r d/dr + f(l)] P_i(r); dr}
!                     
!---------------------------------------------------------------------

      Use RADIAL, L => lro, MX => mro

      Implicit real(8) (A-H,O-Z)
      Integer(4), Intent(in) :: I,J

      DIMENSION Q(NR)

      JJ=I; II=J;    Q(:)=P(:,JJ)*R(:);    LI=L(I); LJ=L(J)

      IL = IABS(LI - LJ)
      IF (IL .NE. 0 .AND. IL .NE. 2) Stop 'GRAD2: LI - LJ <> 0,2'
	  
      A1=(LI+LJ+D2)/((LI+D1)*(LJ+D1))
      A2=((LJ+D5)*(LJ+D1)+(LJ+D1+D5)*(LI+D1))/((LI+D1)*(LJ+D1))
      A = A1 - A2*(LI+LJ+D3)/((LI+LJ+D4)*(LJ+D5))
      FACT=(LJ+D5)/(LI+LJ+D3)
      G=R(1)**2*P(1,I)*P(1,J)*FACT*(1.+A*Z*R(1))

      MM=MIN0(MX(I)+1,MX(J)+1,NR-2)

      K=2
      F1=D5*(P(K+1,II)-P(K-1,II))
      F2=P(K+1,II)-D2*P(K,II)+P(K-1,II)

      G0=Q(K)*R(K)
      G1=D5*(Q(K+1)*R(K+1)-Q(K-1)*R(K-1))
      G2=Q(K+1)*R(K+1)-D2*Q(K)*R(K)+Q(K-1)*R(K-1)

      G=G+D2*F1*G0+(D2*F2*G1+F1*G2)/D3

      DO K=4,MM,2
       F1=D5*(P(K+1,II)-P(K-1,II))
       F2=P(K+1,II)-D2*P(K,II)+P(K-1,II)
       F3=D5*(P(K+1,II)-P(K-2,II))-D2*F1
       F4=P(K+2,II)+P(K-2,II)-D4*(P(K+1,II)+P(K-1,II))+D6*P(K,II)
       G0=Q(K)*R(K)
       G1=D5*(Q(K+1)*R(K+1)-Q(K-1)*R(K-1))
       G2=Q(K+1)*R(K+1)-D2*Q(K)*R(K)+Q(K-1)*R(K-1)
       G3=D5*(Q(K+2)*R(K+2)-Q(K-2)*R(K-2))-D2*G1
       G4=Q(K+2)*R(K+2)+Q(K-2)*R(K-2)-D4*(Q(K+1)*R(K+1) &
               +Q(K-1)*R(K-1))+D6*Q(K)*R(K)
       G=G+D2*F1*G0+(D2*F2*G1+F1*G2)/D3 &
            -(F1*G4-F4*G1+D4*(F2*G3-F3*G2))/90.d0
      End do

      U=QUADR(JJ,II,0)
      G=G-D5*U
      IDELTA=LJ-LI
      Select case(idelta)
       Case(-2); GRAD2=G-(LI-D2)*U
       Case( 0); GRAD2=G+(D1+D5)*U
       Case(+2); GRAD2=G+(LI+D3)*U
      End Select

      END FUNCTION GRAD2

