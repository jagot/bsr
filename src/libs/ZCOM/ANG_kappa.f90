!=======================================================================
      Real(8) Function Cjkj (j1,k,j2)
!=======================================================================
!     Computes the relativistic reduced matrix element
!     (see Eq.15, Grant and Pyper, J.Phys.B9,761,1976)
!
!             k
!     (j1 || C  || j2)  =
!
!                                       (j1   k  j2 )
!         (-1)^(j1+1/2) sqrt([j1][j2])
!                                       (1/2  0 -1/2)
!
!     C(j,0,j) = sqrt([j])
!----------------------------------------------------------------------

      Implicit none
      Integer, intent(in) :: j1,k,j2
      Real(8), external :: Z_3j2

      Cjkj = Z_3j2(j1,1,k+k,0,j2,-1) * &
            sqrt(real((j1+1)*(j2+1))) * (-1)**((j1+1)/2)

      End Function Cjkj


!=======================================================================
      Real(8) Function Ckap (kap1,k,kap2)
!=======================================================================
!     Computes the relativistic reduced matrix element
!                          k
!                (kap1 || C  || kap2)
!
!     see Cjkj as origin
!----------------------------------------------------------------------

      Implicit None
      Integer, intent(in) :: kap1,k,kap2
      Integer, external :: j_kappa
      Real(8), external :: Cjkj
      Integer :: j1,j2

      j1 = j_kappa(kap1)
      j2 = j_kappa(kap1)
      Ckap = Cjkj(j1,k,j2)

      End Function Ckap

!=======================================================================
      Integer Function kappa_lj(l,jj)
!=======================================================================
!     kappa = (l-j)*(2j+1);  jj = 2j
!-----------------------------------------------------------------------
      Integer :: l,jj
      kappa_lj = (2*l-jj)*(jj+1)/2
      End Function kappa_lj


!=======================================================================
      Integer Function l_kappa(kappa)
!=======================================================================
!     l-value for given kappa
!-----------------------------------------------------------------------
      Integer :: kappa

      if(kappa.eq.0) Stop 'l_kappa: kappa=0'
      if(kappa.gt.0) then
       l_kappa =  kappa
      else
       l_kappa = -kappa-1
      end if

      End Function l_kappa


!=======================================================================
      Integer Function j_kappa(kappa)
!=======================================================================
!     j-value for given kappa: j_kappa = 2*j
!-----------------------------------------------------------------------
      Integer :: kappa

      if(kappa.eq.0) Stop 'j_kappa: kappa=0'
      if(kappa.gt.0) then
       j_kappa =  kappa+kappa-1
      else
       j_kappa = -kappa-kappa-1
      end if

      End Function j_kappa


