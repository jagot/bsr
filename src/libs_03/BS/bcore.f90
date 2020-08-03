!======================================================================
      Subroutine BCORE
!======================================================================
!     COMPUTES THE ENERGY OF THE COMMON CLOSED SHELLS  (CORE)
!----------------------------------------------------------------------
      Use spline_atomic, nclosd => kclosd
      Use spline_orbitals, L => LBS

      Implicit none

      Integer :: I,J,K
      Real(8) :: SUMI, SUMJ, CA, CB, CN, TI, TIJ, ck,cl, c1,c2, &
                 TK1,TK2,TK3, UK1,UK2, SN,SN2
      Real(8), external :: ZCB, BHL, RKy, NKy, TKy
      Real(8), parameter :: zero = 0.d0, one = 1.d0, two = 2.d0

      EC = zero

      Do I = 1,NCLOSD

       SUMI = 4*L(I)+2
       TI = RKy(I,I,I,I,0)

       Do K = 2,2*L(I),2
        CA = ZCB(L(i),K,L(i)) * (2*L(i)+1)/(4*L(i)+1)
        TI = TI - CA*RKy(I,I,I,I,k)
       End do

       EC = EC + SUMI*((SUMI-1)*TI - BHL(I,I)) / two

       Do J = 1,I-1
         SUMJ = 4*L(J)+2
         TIJ = RKy(I,J,I,J,0)
         Do K=IABS(L(I)-L(J)),L(I)+L(J),2
           CB = ZCB(L(i),K,L(j))/two
           TIJ = TIJ - CB*RKy(I,J,J,I,k)
         End do
         EC = EC + SUMI*SUMJ*TIJ
       End do

      End do

!----------------------------------------------------------------------
!                                                     o-o  interaction:
      if(ioo.eq.0) Return

      Do I = 1,NCLOSD

        if(L(i).eq.0) Cycle
        SUMI = 4*L(I)+2

        TI = 0.d0
        if(L(I).eq.1) TI = +  8.D0/  5.D0  * NKy(I,I,I,I,0)
        if(L(I).eq.2) TI = + 56.D0/ 21.D0  * NKy(I,I,I,I,0) &
                           + 16.D0/ 21.D0  * NKy(I,I,I,I,2)
        if(L(I).eq.3) TI = + 48.D0/ 13.D0  * NKy(I,I,I,I,0) &
                           + 16.D0/ 13.D0  * NKy(I,I,I,I,2) &
                           + 80.D0/143.D0  * NKy(I,I,I,I,4)
        Do K = 2,2*L(I),2
         CA = ZCB(L(i),K,L(i)) * (2*L(i)+1)/(4*L(i)+1)
         SN = NKy(i,i,i,i,k) 
         SN2 = NKy(i,i,i,i,k-2) 
         c1 =-one*(k+3)/(k+1)/(k+k+3)
         c2 = one*(k-2)/k/(k+k-1)
         TI = TI + CA*(c1*SN-c2*SN2)
        End do

        EC = EC + SUMI*(SUMI-1)*TI/two

        Do J = 1,I-1
         if(L(i).eq.0) Cycle
         SUMJ = 4*L(J)+2
         TIJ = zero
         Do K=IABS(L(I)-L(J)),L(I)+L(J),2
           CB = ZCB(L(i),K,L(j))
           SN = NKy(i,j,j,i,k)
           CN = CB*(L(i)+L(j)+k+2)*(L(i)+L(j)-k)   &
                  *(L(i)-L(j)+k+1)*(L(j)-L(i)+k+1) &
                  /(k+1)/(k+2)
           TIJ = TIJ + CN*SN  

           if(k.eq.0) Cycle
 
           TK1=TKy(i,j,j,i,k+1)-TKy(i,j,j,i,k-1)
           TK2=TKy(j,j,i,i,k+1)-TKy(j,j,i,i,k-1)
           TK3=TKy(i,i,j,j,k+1)-TKy(i,i,i,i,k-1)
           UK1 = TK1 + TK2
           UK2 = TK1 + TK3

           ck = k*(k+1)
           cl = L(i)*(L(i)+1)-L(j)*(L(j)+1)
           c1 = cl - ck
           c2 =-cl - ck 
           TIJ = TIJ + CB*(ck*TK1 + (c1*UK1 + c2*UK2)/two)

           SN2 = NKy(i,j,j,i,k-2)
           cn = CB*c1*c2/two
           c1 =-cn*(k+3)/(k+1)/(k+k+3)
           c2 = cn*(k-2)/k/(k+k-1)
           TIJ = TIJ + c1*SN + c2*SN2 
         End do

         EC = EC + SUMI*SUMJ*TIJ

       End do

      End do

      End Subroutine BCORE
