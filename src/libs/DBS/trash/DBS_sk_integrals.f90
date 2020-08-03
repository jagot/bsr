!====================================================================
      Module DBS_sk_integrals
!====================================================================
!     contains the B-spline representation of two-electron integrals
!     rkb(i,j;i',j') in the non-symmetric column storage mode:
!                rkb(1:ns, 1:ns, 1:2*ks-1, 1:2*ks-1)
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

      Real(8) :: memory_DBS_sk_integrals = 0.d0

      End Module DBS_sk_integrals


!====================================================================
      Subroutine alloc_DBS_sk_integrals(ns,ks,kmin,kmax,ktype)
!====================================================================
! ... Allocates space for spline integrals:
!     1. ppqq
!     2. pqqp
!--------------------------------------------------------------------
      Use DBS_sk_integrals

      Implicit none
      Integer, intent(in) :: ns,ks,kmin,kmax,ktype

      if(allocated(irka)) Deallocate (irka,rka)
      Allocate(irka(kmin:kmax,ktype))
      Allocate(rka(ns,ns,2*ks-1,2*ks-1,kmin:kmax,ktype))
      kra_min = kmin
      kra_max = kmax
      irka = -100
      ntype = ktype

      if (associated(rkb)) Nullify(rkb)
      itype='bbbb';  krk=-100

      memory_DBS_sk_integrals = ((kmax-kmin+1)*ntype*4 + &
       ns*ns*ks*ks*(kmax-kmin+1)*ktype*8)/(1024d0*1024d0)

      Call alloc_DBS_moments

      End Subroutine alloc_DBS_sk_integrals


!====================================================================
      Subroutine dealloc_DBS_sk_integrals
!====================================================================
! ... deallocates arrays in module "DBS_sk_integrals"
!--------------------------------------------------------------------
      Use DBS_sk_integrals

      if (associated(rkb)) nullify(rkb)
      itype='aaaa'; krk = -100
      if (allocated(irka)) Deallocate(irka,rka)
      kra_min = -100
      kra_max = -100
      ntype = 0
      memory_DBS_sk_integrals = 0.d0

      End Subroutine dealloc_DBS_sk_integrals



