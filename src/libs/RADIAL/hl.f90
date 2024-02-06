!------------------------------------------------------------------
      Real(8) FUNCTION HL(I,J)
!------------------------------------------------------------------
!
!     returns the value of
!
!     <i|L|j> =  < P_i(r) | d^2/dr^2 + 2Z/r - l(l+1)/r^2 | P_j(r) >
!
!     using a special formula to preserve symmetry
!------------------------------------------------------------------

      USE RADIAL, L => lro

      IMPLICIT REAL(8) (A-H,O-Z)

      INTEGER(4), INTENT(in) :: I,J

      IF (IABS(L(I)-L(J)).NE.0) Stop ' HL: L_i <> L_j '

      LI = L(I)
      C  = 2*LI + 1
      A1 = -D2/(C*(LI+1))
      A2 = A1/((C+D2)*(LI+1))
      A3 = A2/((LI+2)*(LI+1))
      ZR = Z*R(1)
      HL = H*C*P(1,I)*P(1,J)*(D1+ZR*(A1+ZR*(A2+ZR*A3)))

      K = 2
      C = D4/D3
      DI1 = P(K+1,I) - P(K-1,I)
      DI2 = P(K+1,I) - D2*P(K,I) + P(K-1,I)
      DJ1 = P(K+1,J) - P(K-1,J)
      DJ2 = P(K+1,J) - D2*P(K,J) + P(K-1,J)
      HL = HL + DI1*DJ1 + C*DI2*DJ2

      MM = MIN0(mro(I)+3,mro(J)+3,NR-3)
      DO K = 4,MM,2
      DI1 = P(K+1,I) - P(K-1,I)
      DI2 = P(K+1,I) - D2*P(K,I) + P(K-1,I)
      DI4 = P(K+2,I) - D4*(P(K+1,I)+P(K-1,I)) + D6*P(K,I) +P(K-2,I)
      DI3 = P(K+2,I) - P(K-2,I) - D2*DI1
      DI5 = P(K+3,I)-P(K-3,I) - D4*(P(K+2,I)-P(K-2,I)) &
          + 5.D0*(P(K+1,I)-P(K-1,I))
      DI6 = P(K+3,I)+P(K-3,I) - D6*(P(K+2,I)+P(K-2,I)) &
          + 15.D0*(P(K+1,I)+P(K-1,I)) - 20.D0*P(K,I)
      DJ1 = P(K+1,J) - P(K-1,J)
      DJ2 = P(K+1,J) - D2*P(K,J) + P(K-1,J)
      DJ4 = P(K+2,J) - D4*(P(K+1,J)+P(K-1,J)) + D6*P(K,J) +P(K-2,J)
      DJ3 = P(K+2,J) - P(K-2,J) - D2*DJ1
      DJ5 = P(K+3,J)-P(K-3,J) - D4*(P(K+2,J)-P(K-2,J)) &
            + 5.D0*(P(K+1,J)-P(K-1,J))
      DJ6 = P(K+3,J)+P(K-3,J) - D6*(P(K+2,J)+P(K-2,J)) &
            + 15.D0*(P(K+1,J)+P(K-1,J)) - 20.D0*P(K,J)

      HL = HL + DI1*DJ1 + C*DI2*DJ2 + (DI3*DJ3 + DI2*DJ4+DI4*DJ2)/45.D0 &
         -(DI3*DJ5+DI5*DJ3)/252.D0 - (DI2*DJ6+DI6*DJ2-1.1*DI4*DJ4)/378.D0

      End do

      TZ = Z + Z
      C = (LI + D5)**2
      HL2 = D5*(TZ*R(1) - C)*P(1,I)*P(1,J)
      DO K = 2,MM,2
       HL2 = HL2 + D2*(TZ*R(K  ) - C)*P(K  ,I)*P(K,  J) &
                 +    (TZ*R(K+1) - C)*P(K+1,I)*P(K+1,J)
      END DO

      HL = -HL/(D2*H) + HL2*H1

      IF(rel) HL = HL - D2*RLSHFT(I,J)

      END FUNCTION HL

