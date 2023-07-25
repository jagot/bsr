!=======================================================================
      Integer Function kappa_lj(l,jj)
!=======================================================================
!     l, j  ->  kappa = (l-j)*(2j+1);  jj = 2j
!-----------------------------------------------------------------------
      Integer :: l,jj
      kappa_lj = (2*l-jj)*(jj+1)/2
      End Function kappa_lj

!=======================================================================
      Integer Function l_kappa(kappa)
!=======================================================================
!     kappa ->  l  
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
!     kappa ->  j   (= 2*j)
!-----------------------------------------------------------------------
      Integer :: kappa
      if(kappa.eq.0) Stop 'j_kappa: kappa=0'
      if(kappa.gt.0) then
       j_kappa =  kappa+kappa-1
      else
       j_kappa = -kappa-kappa-1
      end if 
      End Function j_kappa

