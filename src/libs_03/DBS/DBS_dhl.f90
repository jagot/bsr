!====================================================================
      Module DBS_dhl_pq
!====================================================================
!     contains the spline representation of Dirac hamiltonian
!      <P|V_nucl|P>     c<P|-d/dr+k/r|Q>
!     c<Q|d/dr+k/r|P>   <Q|V_nucl-2c^2|Q> 
!--------------------------------------------------------------------
      Implicit none

      Integer :: k_dhl = 0               ! k-value for HD-operator
      Real(8), allocatable :: dhl(:,:)

! ... Bloch operator:

      Integer :: mbloch  =   1    !  flag to include  
      Real(8) :: RA = 0.d0        !  boder radius (=tmax)
      Real(8) :: RB = 0.d0        !  Bloch b-parameter
      Real(8) :: RN = 0.d0        !  Q(a)/P(a)
      Real(8) :: pnu = 1.d0       !  Q(a)/P(a)

      End Module DBS_dhl_pq

!======================================================================
      Subroutine alloc_dhl_pq(m)
!======================================================================
!     allocate or deallocate dhl array
!----------------------------------------------------------------------
      Use DBS_dhl_pq

      Integer :: m
      if(m.eq.0) then
       if(allocated(dhl)) Deallocate(dhl)
       k_dhl=0
      else
       if(allocated(dhl)) Deallocate(dhl)
       Allocate(dhl(m,m))
       k_dhl=0; dhl = 0.d0
      end if

      End Subroutine alloc_dhl_pq

!======================================================================
      Subroutine MAT_dhl_pq(k)
!======================================================================
!     Builds full matrix for the coulomb problem with kappa=k
!     from the symmetric forms for elementary operators
!----------------------------------------------------------------------
      USE zconst, only: c_au
      USE DBS_grid
      USE DBS_gauss
      USE DBS_dhl_pq

      Implicit none
      Real(8) :: sa,sb
      Integer :: k,i,j

      if(k_dhl.eq.k) Return
      if(.not.allocated(dhl)) Call Alloc_dhl_pq(ms)

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

      k_dhl = k

! ... symmetrize dhl by Bloch operator: 

      if(mbloch.eq.1) then
       RA = tmax
       RN = (RB+k)/(2.d0*RA*c_au)
       i = nsp; j = ns + nsq
       dhl(i,i) = dhl(i,i) - pnu*c_au*RN
       dhl(i,j) = dhl(i,j) + pnu*c_au  
       dhl(j,i) = dhl(j,i) + (pnu-1.d0)*c_au
       dhl(j,j) = dhl(j,j) - (pnu-1.d0)*c_au/RN 
      end if

      End Subroutine MAT_dhl_pq


!======================================================================
      Real(8) Function Vp_dhl(k1,v1,k2,v2)
!======================================================================
!     integral   <v1| dhl |v2>
!----------------------------------------------------------------------
      USE DBS_grid,   only: ms
      USE DBS_dhl_pq

      Implicit none
      Integer, Intent(in) :: k1,k2
      Real(8), Intent(in) :: v1(ms),v2(ms)
      Real(8) :: v(ms)

      Vp_dhl = 0.d0; if(k1.ne.k2) Return
      Call Mat_dhl_pq(k1)
      v = MATMUL(dhl,v2)
      Vp_dhl = SUM(v1*v)

      End Function Vp_dhl



