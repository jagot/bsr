!=======================================================================
      REAL(8) FUNCTION RLSHFT (i1,i2)
!=======================================================================
!
!     relativistic shift:  mass correction and one-body Darwin terms
!
!-----------------------------------------------------------------------

      USE RADIAL, L => lro, MX => mro

      IMPLICIT REAL(8) (A-H,O-Z)

!     FORM  DD - L(L+1)/RR |P(I)>

      FL = L(I1)
      C  = (FL+D5)**2
      LL = L(I1) + 1
      L2 = 2*L(I1) + 1
      L3 = 2*L(I1) + 3
      Z2 = Z*Z
      HH = 180.D0*H*H
      MM = MAX0(MX(I1),MX(I2))
      Do J = 2,MM
       YK(J) = -D1/RR(J)
      End do

!     FORM THE INTEGRAND

      I = I1
      A1 = D0
      B1 = D0

      Do KK = 1,2
      B2 = B1
      YY = (P(3,I)+P(1,I) - D2*P(2,I))/(H*H) - C*P(2,I)
      YK(2) = YY*YK(2)
      YK(3) = YK(3)*((-(P(5,I)+P(1,I)) + D16*(P(4,I)+P(2,I)) &
              -D30*P(3,I))/(D12*H*H) - C*P(3,I))
       MM = MX(I) - 3
       Do  K =  4,MM
        YY = D2*(P(K+3,I)+P(K-3,I))
        YY = YY - 27.D0*(P(K+2,I)+P(K-2,I))
        YY = YY +  270.D0*(P(K+1,I)+P(K-1,I)) - 490.D0*P(K,I)
        YY = YY/HH - C*P(K,I)
        YK(K) = YY*YK(K)
        IF (K .EQ. 4) B1 = (YY/(D2*Z*P(4,I)*R(4)) + D1)/R(4)
       End do

      MM = MM + 1
      YK(MM) = YK(MM)*((-(P(MM+2,I)+P(MM-2,I))+D16*(P(MM+1,I)+P(MM-1,I)) &
              -D30*P(MM,I))/(D12*H*H) - C*P(MM,I))
      MM = MM + 1
      YK(MM) = YK(MM)*((P(MM+1,I) + P(MM-1,I) - D2*P(MM,I))/(H*H) &
             - C*P(MM,I))
      A2 = A1
      A1 = (P(1,I)/(AZ(I)*R(1)**L(I)*R2(1)) - D1 +Z*R(1)/LL)/RR(1)
      I = I2
      End do

! ... DETERMINE CONTRIBUTION FROM NEAR THE NUCLEUS

      A = (Z/LL - L2*(B1 + B2)/D2)/LL
      B = (L2*B1*B2 - D2*(A1 + A2) + (Z/LL**2)*(D2*Z*(D1 + D1/LL) &
           - L2*(B1 + B2)))/L3
      RELSH = -P(4,I1)*P(4,I2)*(D1 + A*R(4) + B*RR(4))*D4*Z2/L2

      RELSH = RELSH/H1 - D5*YK(4)

      RELSH2 = D0

! ... INTEGRATE

      MM = MAX0(MX(I1),MX(I2))
      Do J = 5,MM,2
       RELSH2 = RELSH2 + YK(J)
       RELSH  = RELSH  + YK(J-1)
      End do

      RELSH  = (RELSH + D2*RELSH2)*H1

! ... add one-body Darwin term:

      IF (L(i1).eq.0) RELSH = RELSH + Z*AZ(i1)*AZ(i2)

      RLSHFT = RELSH*D5*FINE

      END FUNCTION RLSHFT
