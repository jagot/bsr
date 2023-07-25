!======================================================================
      SUBROUTINE GET_ECORE(ecore)
!======================================================================
!     COMPUTES THE ENERGY OF THE COMMON CLOSED SHELLS ( = nclosd )
!----------------------------------------------------------------------

      USE RADIAL

      IMPLICIT NONE
      Real(8), Intent(out) :: ecore
      Integer(4) :: k,i,j, k1,k2, m1,m2
      Real(8) :: ca,cb,cn, TI,TJ
      Real(8), External :: RK, NK, HL, ZCB, RK_oo

      ECORE = 0.D0

      DO I = 1,kCLOSD

! ... intra-shell interaction ...

       TI = RK(I,I,I,I,0)
       DO K = 2,2*lro(I),2
        CA = - ZCB(lro(i),K,lro(i)) * (2*lro(i)+1)/(4*lro(i)+1)
        TI = TI + CA*RK(I,I,I,I,k)
       END DO

       if(oo) then                                ! o-o contribution
        if(lro(I).eq.1) then
          TI = TI +  8.D0/  5.D0  * NK(I,I,I,I,0)
        elseif(lro(I).eq.2) then
          TI = TI + 56.D0/ 21.D0  * NK(I,I,I,I,0)
          TI = TI + 16.D0/ 21.D0  * NK(I,I,I,I,2)
        elseif(lro(I).eq.3) then
          TI = TI + 48.D0/ 13.D0  * NK(I,I,I,I,0)
          TI = TI + 16.D0/ 13.D0  * NK(I,I,I,I,2)
          TI = TI + 80.D0/143.D0  * NK(I,I,I,I,4)
        end if
       end if

       ECORE = ECORE + (4*lro(I)+2)*( (4*lro(I)+1)*TI - Hl(I,I) ) / D2

!  ... inter-shell interaction ...

       DO J = 1,I-1

         TJ = RK(I,J,I,J,0)
         DO K = IABS(lro(I)-lro(J)), lro(I)+lro(J), 2
           CB = - ZCB(lro(i),K,lro(j))/2
           TJ = TJ + CB*RK(I,J,J,I,k)

           if(oo) then                        ! o-o contribution
            TJ = TJ + CB*RK_oo(I,J,J,I,k)
            if(lro(i).gt.0.and.lro(j).gt.0) then
             k1 = k + 1; k2 = k + 2
             m1 = lro(i) + lro(j) + 1;  m2 = lro(i) - lro(j)
             CN = CB * (m1*m1-k1*k1)*(k1*k1-m2*m2)/(k1*k2)
             if(abs(CN).gt.0.d0) then
               TJ = TJ + CN*(NK(I,J,J,I,k) + NK(J,I,I,J,k))/D2
             end if
            end if
           end if

         END DO

         ECORE = ECORE + (4*lro(I)+2)*(4*lro(J)+2)*TJ
       END DO
      END DO

      END SUBROUTINE GET_ECORE