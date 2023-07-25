!======================================================================
      Subroutine TVM(k,p,q,TA,VA,MA)
!======================================================================
!     Calculate the expectation values for 
!     TA - kinatic energy
!     VA - potential energy
!     MA - mass energy
!     for orbital (p,q) with kappa value "k"
!----------------------------------------------------------------------
      USE zconst
      USE DBS_nuclear, z => atomic_number
      USE DBS_grid
      USE DBS_gauss

      Implicit none
      Integer, intent(in) :: k
      Real(8), intent(in) :: p(ns),q(ns)
      Real(8), intent(out) :: TA,VA,MA
      Real(8) :: sa,sb,vv(ns),a(ns,ns)

! ... potential average:

      vv=MATMUL(fpb_nucl,p)
      VA = dot_product(p,vv)
      vv=MATMUL(fqb_nucl,q)
      VA = VA + dot_product(q,vv)

! ... mass average:

      sb = -2.d0 * c_au * c_au
      vv=MATMUL(fqbs,q)
      MA = sb*dot_product(q,vv)

! ... kinetic average:

! ... matrix of (2:1) = c <B | D+ | B> 

      sa = c_au; sb = dble(k)*c_au
      a(1:ns,1:ns) = sa*fqpbsd + sb*fqpbs
      vv=MATMUL(a,p)
      TA = dot_product(q,vv)

! ... matrix of (1:2) = c <B | D- | B> 

      sa = -c_au; sb = dble(k)*c_au
      a(1:ns,1:ns) = sa*fpqbsd + sb*fpqbs
      vv=MATMUL(a,q)
      TA = TA + dot_product(p,vv)

      End Subroutine TVM

