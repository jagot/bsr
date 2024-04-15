!====================================================================
      Real(8) Function ZCKLM (L1,M1,L2,M2,K)
!====================================================================
!     Condon-Shortly coefficients: ???
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: L1,M1,L2,M2,K
      Real(8), external :: CLEBSH
      Integer :: I,I1,I2, J,J1,J2, IM

      ZCKLM=0.0
      if(K.GT.L1+L2.OR.K.LT.IABS(L1-L2)) Return
      I=L1+L2+K
      if(I.NE.I/2*2) Return

      J = K+K+1
      J1= L1+L1+1
      J2= L2+L2+1
      I1=-M1-M1+1
      I2= M2+M2+1
      IM= I1+I2-1

      ZCKLM = SQRT(DBLE(J1*J2))/J*(-1)**M2  *  &
              CLEBSH(J1,1,J2,1,J,1)*CLEBSH(J1,I1,J2,I2,J,IM)

      End Function ZCKLM
