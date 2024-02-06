!=======================================================================
      Real(8) function Cjkj (j1,k,j2)
!=======================================================================
!     Computes the relativistic reduced matrix element
!     (see Eq.15, Grant and Pyper, J.Phys.B9,761,1976)
!
!     (j1 || C^k  || j2)  =
!                                       (j1   k  j2 )
!         (-1)^(j1+1/2) sqrt([j1][j2])
!                                       (1/2  0 -1/2)
!     C(j,0,j) = sqrt([j])
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j1,k,j2
      Real(8), external :: Z_3j2
      Cjkj = Z_3j2(j1,1,k+k,0,j2,-1) * &
            sqrt(real((j1+1)*(j2+1))) * (-1)**((j1+1)/2)
      End function Cjkj

!=======================================================================
      Real(8) function Ckap (kap1,k,kap2)
!=======================================================================
!     Computes the relativistic reduced matrix element
!
!                (kap1 || C^k  || kap2)
!
!     see Cjkj as origin
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: kap1,k,kap2
      Integer, external :: j_kappa
      Real(8), external :: Cjkj
      Integer :: j1,j2

      j1 = j_kappa(kap1);  j2 = j_kappa(kap2);  Ckap = Cjkj(j1,k,j2)

      End function Ckap
