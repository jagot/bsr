!======================================================================
      Real(8) FUNCTION  ZETA(I1,I2)
!======================================================================
!
!     COMPUTES THE NUCLEAR SPIN-ORBIT PARAMETER
!
!     alfa^2/4   Int(r,0,inf) [ P_i1(r)  z/r^3  P_i2(r) ]
!
!     AND THE CORRECTIONS FOR THE COMMON CLOSED SHELLS
!
!----------------------------------------------------------------------

      USE RADIAL, L => lro

      IMPLICIT NONE
      Integer(4), Intent(in) :: I1,I2
      Integer(4) :: LA,LB,IN,IM,IV,IW,I,K,K1,KK,M1,M2
      REAL(8) :: CA,CB,CN,CM,CV
      Real(8), External :: NK, VK, QUADR, ZCB

      if(L(i1).ne.L(i2)) Stop ' ZZETA: L(i1) <> L(i2)'
      if(L(i1).eq.0    ) Stop ' ZZETA: L = 0 '

      ZETA = FINE*Z*QUADR(I1,I2,-3)

      LB = L(I1)
      Do I = 1,KCLOSD
         LA = L(I)
         ca = 4*LA+2
         ZETA = ZETA - ca*NK(I1, I, I2, I, 0)

         Do k=iabs(LB-LA),LB+LA,2
          cb=3*ca*ZCB(LA,k,LB)/(8*LB*(LB+1))
          CN = (NK(I1,I,I,I2,K) + NK(I,I1,I2,I,K))/(k+1)
          CM = D0
          CV = D0
          if(k.gt.0) then
           CM = (NK(I1,I,I,I2,K-2) + NK(I,I1,I2,I,K-2))/k
           CV =  2   *(VK(I1,I,I,I2,K-1) - VK(I,I1,I2,I,K-1)) + &
                 k   *(NK(I1,I,I,I2,K  ) - NK(I,I1,I2,I,K  )) + &
                (k+1)*(NK(I1,I,I,I2,K-2) - NK(I,I1,I2,I,K-2))
          end if

          IW = LB*(LB+1)-LA*(LA+1);  IV = IW + k*(k+1)
          M1 = LA+LB+1; M1=M1*M1
          M2 = LA-LB;   M2=M2*M2
          K1 = (K+1)*(K+1); KK=K*K
          IN = (M1-K1)*(K1-M2)
          IM = (M1-KK)*(KK-M2)

          ZETA = ZETA + cb*(IV*(CV + IW*(CN-CM)) + IN*CN - IM*CM)

         End do    !  over k
       End do      !  over core orbitals

      END FUNCTION ZETA
      