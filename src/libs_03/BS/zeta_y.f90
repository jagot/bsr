!======================================================================
      Real(8) FUNCTION ZETA_y (I1,I2)
!======================================================================
!
!   COMPUTES THE NUCLEAR SPIN-ORBIT PARAMETER AND THE
!   CORRECTIONS FOR THE COMMON CLOSED SHELLS
!
!   Calls:  QUADR, NKy, VKy  - integrals
!           ZCB - 3j-symbol (to determine angular part of the intershell
!                            interaction) 
!----------------------------------------------------------------------
      USE spline_orbitals, p => pbs, L => lbs
      USE spline_atomic, nclosd => kclosd

      IMPLICIT NONE

      INTEGER, INTENT(in) :: i1,i2

      REAL(KIND=8), EXTERNAL :: QUADR, NKy, VKy, ZCB

      INTEGER :: i, LA,LB, k,kmin,kmax, k1,k2, m1,m2, iv,iw, ink,in2
      REAL(KIND=8) :: ZETA, ZZ, QA, CB
      REAL(KIND=8) :: DNK,ENK,SNK,DVK,EVK,SVK,SN2,DN2,EN2

      if(L(i1).ne.L(i2)) Stop ' ZETA: L(i1) <> L(i2)'
      if(L(i1).eq.0) Stop ' ZETA: L = 0 '

      ZETA = fine*Z*QUADR(I1,I2,-3)

      LB = L(I1)
      DO I = 1,NCLOSD
        LA = L(I)
        QA = 4*LA+2
        ZETA = ZETA - QA*NKy(I1, I, I2, I, 0)

        kmin=iabs(LB-LA)
        kmax=iabs(LB+LA)
        Do k=kmin,kmax,2

          cb = 3*QA*ZCB(LA,k,LB)/(8*LB*(LB+1))

          k1 = k + 1
          DNK = NKy(I1, I, I, I2, K)
          ENK = NKy(I, I1, I2, I, K)
          SNK = (DNK + ENK)/(k+1)
          SVK = 0.d0
          SN2 = 0.d0

          if(k.gt.0) then
           DN2 = NKy(I1, I, I, I2, K-2)
           EN2 = NKy(I, I1, I2, I, K-2)
           SN2 = (DN2 + EN2)/k
           DVK = 2*VKy(I1,I,I,I2,K-1) - k*DNK + (k+1)*EN2
           EVK = 2*VKy(I,I1,I2,I,K-1) - k*ENK + (k+1)*DN2
           SVK = DVK - EVK
          end if

          IW=LB*(LB+1)-LA*(LA+1)
          IV=IW+k*k1
          M1=(LA+LB+1)
          m1=m1*m1
          m2=m2*m2
          M2=(LA-LB)
          K1=(K+1)*(K+1)
          K2=K*K
          INK=(M1-K1)*(K1-M2)
          IN2=(M1-K2)*(K2-M2)

          ZZ = IV*(SVK + IW*(SNK-SN2)) + INK*SNK - IN2*SN2
          ZETA = ZETA + ZZ*cb

        End do
      End do

      ZETA_y = ZETA

      End Function ZETA_y

