!====================================================================
      Module DBS_integrals
!====================================================================
!     contains the B-spline representation of two-electron integrals
!     rkb(i,j;i',j') in the symmetric column storage mode:
!                rkb(1:ns, 1:ns, 1:ks, 1:ks)
!     itype - character which indicates the type of integral
!     krk   - multipole index for the integral
!--------------------------------------------------------------------
      Implicit none

      Integer :: krk = -100
      Real(8), pointer :: rkb(:,:,:,:)
      Character(4) :: itype='aaaa'

      Integer :: kra_min = -100
      Integer :: kra_max = -100
      Integer :: ntype = 0
      Integer, allocatable :: irka(:,:)
      Real(8), allocatable, target :: rka(:,:,:,:,:,:)

      Real(8) :: memory_DBS_integrals = 0.d0

      End Module DBS_integrals


!====================================================================
      Subroutine alloc_DBS_integrals(ns,ks,kmin,kmax,ktype)
!====================================================================
! ... Allocates space for spline integrals:
!     1.  pppp  
!     2.  qqqq  
!     3.  pqpq
!     4.  qpqp
!--------------------------------------------------------------------
      Use DBS_integrals

      Implicit none
      Integer, intent(in) :: ns,ks,kmin,kmax,ktype

      if(allocated(irka)) Deallocate (irka,rka)
      Allocate(irka(kmin:kmax,ktype))
      Allocate(rka(ns,ns,ks,ks,kmin:kmax,ktype))
      kra_min = kmin
      kra_max = kmax
      irka = -100
      ntype = ktype

      if (associated(rkb)) Nullify(rkb)
      itype='bbbb';  krk=-100

      memory_DBS_integrals = ((kmax-kmin+1)*ntype*4 + &
       ns*ns*ks*ks*(kmax-kmin+1)*ktype*8)/(1024d0*1024d0)

      Call alloc_DBS_moments

      End Subroutine alloc_DBS_integrals


!====================================================================
      Subroutine dealloc_DBS_integrals
!====================================================================
! ... deAllocates arrays in module "spline_integrals"
!--------------------------------------------------------------------
      Use DBS_integrals

      if (associated(rkb)) nullify(rkb)
      itype='aaaa'; krk = -100
      if (allocated(irka)) Deallocate(irka,rka)
      kra_min = -100
      kra_max = -100
      ntype = 0
      memory_DBS_integrals = 0.d0

      End Subroutine dealloc_DBS_integrals



