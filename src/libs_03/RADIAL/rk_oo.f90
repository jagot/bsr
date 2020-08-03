!====================================================================
      REAL(8) FUNCTION RK_oo (i1,i2,i3,i4,k)
!====================================================================
!
!     Correction to RK integral due to o-o interaction 
!
!==================================================================

      USE RADIAL, L => lro

      IMPLICIT NONE
      Integer(4), Intent(in) :: i1,i2,i3,i4,k
      REAL(8) :: T,U1,U2,M1,M2, C,C1,C2
      REAL(8), EXTERNAL :: TK,UK,NK

      RK_oo = D0
      if(k.eq.0) Return
      if(i1.eq.i2.and.i1.eq.i3.and.i1.eq.i4) Return

      T  = TK(i1,i2,i3,i4,k+1) - TK(i1,i2,i3,i4,k-1)
      U1 = UK(i1,i2,i3,i4,k+1) - UK(i1,i2,i3,i4,k-1)
      U2 = UK(i2,i1,i4,i3,k+1) - UK(i2,i1,i4,i3,k-1)

      M1 = NK(i1,i2,i3,i4,k-2) + NK(i2,i1,i4,i3,k-2)
      M2 = NK(i1,i2,i3,i4,k  ) + NK(i2,i1,i4,i3,k  )

      C = 2*k*(k+1)
      C1 = L(i1)*(L(i1)+1) - L(i3)*(L(i3)+1) - k*(k+1)
      C2 = L(i2)*(L(i2)+1) - L(i4)*(L(i4)+1) - k*(k+1)
      RK_oo = - C*T - C1*U1 - C2*U2

      C  = C1*C2/D2
      C1 = C*(k-2)/k/(k+k-1)
      C2 = C*(k+3)/(k+1)/(k+k+3)/D2                     
      RK_oo = RK_oo - C1*M1 + C2*M2

      END function RK_oo

