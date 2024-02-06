!====================================================================
      Real(8) Function T_LS_jk(l1,l2,s1,s2,L,S,js,jk,J)
!====================================================================
!     the recoupling coefficient from LS- to jK-coupling:
!
!     < (l1,l2)L,(s1,s2)S;J | (((l1,s1)js,l)K,s2)J >
!
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l1,l2,s1,s2,L,S,js,jk,J
      Real(8), external :: Z_6j

      T_LS_jk = Z_6j(L,s1,jk,s2,J,S) * Z_6j(l2,l1,L,s1,jk,js) * &
                sqrt(1.d0*L*S*js*jk)* (-1)**((s2+J-l2-js)/2)

      End Function T_LS_jk


!====================================================================
      Real(8) Function T_LS_jj(l1,l2,s1,s2,L,S,j1,j2,J)
!====================================================================
!     recoupling coefficient from LS- to jj-coupling:
!
!     < (l1,l2)L,(s1,s2)S;J | ((l1,s1)j1,(l2,s2)j2;J >
!
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l1,l2,s1,s2,L,S,j1,j2,J
      Integer :: k,k_min,k_max
      Real(8), external :: Z_6j

      T_LS_jj = 0.d0
      k_min = max(iabs(l2-j1),iabs(s2-J),iabs(L-s1))+1
      k_max = max(iabs(l2+j1),iabs(s2+J),iabs(L+s1))-1
      if(k_min.gt.k_max) Return
      Do k = k_min,k_max,2
       T_LS_jj = T_LS_jj + Z_6j(l2,s2,j2,J,j1,k) * &
                           Z_6j(L,S,J,s2,k,s1)   * &
                           Z_6j(l1,s1,j1,K,l2,L) * &
                           k * (-1)**(k-1)
      End do
      T_LS_jj = T_LS_jj * sqrt(1.d0*L*S*j1*j2)

      End Function T_LS_jj


!====================================================================
      Real(8) Function T_jj_jk(l2,s2,j1,j2,K,J)
!====================================================================
!     the recoupling coefficient from jK- to jj-coupling:
!
!     < (((l1,s1)j1,l2)K,s2)J | (((l1,s1)j1,(l2,s2)j2,J >
!
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l2,s2,j1,j2,K,J
      Real(8), external ::  Z_6j2

      T_jj_jk = Z_6j2(j1,l2,K,s2,J,j2) * sqrt(1.d0*(K+1)*(j2+1)) &
                *(-1)**((j1+l2+s2+J)/2)

      End Function T_jj_jk
