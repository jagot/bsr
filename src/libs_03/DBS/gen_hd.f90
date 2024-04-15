!======================================================================
      Subroutine Gen_hd(k,dhl)
!======================================================================
!     Builds full matrix for the coulomb problem with kappa=k
!     from the symmetric forms for elementary operators
!----------------------------------------------------------------------
      Use zconst, only: c_au
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: k
      Real(8) :: dhl(ms,ms),sa,sb

! ... matrix of (1:1) = <B|V|B>
      dhl(1:ns,1:ns) = fpb_nucl
! ... matrix of (2:2) = <B|V|B> - 2*c^2*<B|B>
      sb = -2.d0 * c_au * c_au
      dhl(ns+1:ms,ns+1:ms) = fqb_nucl + sb*fqbs
! ... matrix of (2:1) = c <B | D+ | B>
      sa = c_au; sb = dble(k)*c_au
      dhl(ns+1:ms,1:ns) = sa*fqpbsd + sb*fqpbs
! ... matrix of (1:2) = c <B | D- | B>
      sa = -c_au; sb = dble(k)*c_au
      dhl(1:ns,ns+1:ms) = sa*fpqbsd + sb*fpqbs

      End Subroutine Gen_hd
