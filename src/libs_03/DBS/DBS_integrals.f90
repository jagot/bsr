!====================================================================
      Module DBS_integrals
!====================================================================
!     contains the B-spline representation of two-electron integrals
!     rkb(i,j;i',j') in the column storage mode:
!             rkb(1:ns, 1:ns, 1:ks, 1:ks)  - symmetric
!             rkb(1:ns, 1:ns, 1:2*ks1, 1:2*ks+1)  - non-symmetric
!     itype - character which indicates the type of integral
!     krk   - multipole index for the integral
!--------------------------------------------------------------------
      Implicit none

      Integer :: krk = -100
      Real(8), pointer :: rkb(:,:,:,:)
      Character(4) :: itype='aaaa'

      ! storage for Rk-integrals if any

      Integer :: kra_min = -100
      Integer :: kra_max = -100
      Integer :: ntype = 0
      Integer, allocatable :: irka(:,:)
      Real(8), allocatable, target :: rka(:,:,:,:,:,:)
      Real(8), allocatable :: rkb_tmp(:,:,:,:)

      Real(8), allocatable, target :: rka1(:,:,:,:)
      Integer :: ntype1 = 0
      Integer :: krk1 = -100
      Character(4) :: itype1='aaaa'

      ! Storage for Rk-integrals partially contracted with density
      ! matrices. The integers after "d" indicate which indices of the
      ! Rk-integrals the density matrix is contracted with, for use
      ! with this particular storage.
      Real(8), allocatable :: d24_rk(:,:,:)
      Real(8), allocatable :: d23_rk(:,:,:)

      ! storage for Sk-integrals if any

      Integer :: krs_min = -100
      Integer :: krs_max = -100
      Integer :: stype = 0
      Integer, allocatable :: isk(:,:)
      Real(8), allocatable, target :: rks(:,:,:,:,:,:)

      Real(8), allocatable, target :: rks1(:,:,:,:)
      Integer :: jtype1 = 0
      Integer :: krs1 = -100
      Character(4) :: stype1='aaaa'

      Real(8) :: memory_DBS_integrals = 0.d0
      Real(8) :: memory_Rk_integrals  = 0.d0
      Real(8) :: memory_Rk_integral   = 0.d0
      Real(8) :: memory_Sk_integrals  = 0.d0
      Real(8) :: memory_Sk_integral   = 0.d0

      End Module DBS_integrals


!====================================================================
      Subroutine alloc_Rk_integrals(ns,ks,kmin,kmax,ktype)
!====================================================================
!     Allocates/deallocate space for B-spline Rk-integrals
!     for different types and multipole indexes.
!     We expecting four types:   1. pppp  2. qqqq  3. pqpq  4. qpqp
!--------------------------------------------------------------------
      Use DBS_integrals

      Implicit none
      Integer, intent(in) :: ns,ks,kmin,kmax,ktype

      if(associated(rkb)) Nullify(rkb)
      krk=-100
      itype='aaaa'

      if(allocated(irka)) then
       Deallocate (irka,rka)
       ntype = 0
       memory_Rk_integrals = 0.d0
      end if

      if(ktype.gt.0) then
       Allocate(irka(kmin:kmax,ktype))
       Allocate(rka(ns,ns,ks,ks,kmin:kmax,ktype))
       Allocate(rkb_tmp(1:ns,1:ns,1:ks,1:ks))
       Allocate(d24_rk(ns,ks,ks))
       Allocate(d23_rk(ns,ns,ks))
       kra_min = kmin
       kra_max = kmax
       irka = -100
       ntype = ktype
       memory_Rk_integrals = ((kmax-kmin+1)*ntype*4 + &
        ns*ns*ks*ks*(kmax-kmin+1)*ktype*8 + &
        ns*ns*ks*ks*8 + &
        ns*ks*ks*8 + &
        ns*ns*ks*8)/(1024d0*1024d0)
       Call alloc_DBS_moments
      end if

      memory_DBS_integrals =    memory_Rk_integrals  &
                              + memory_Rk_integral   &
                              + memory_Sk_integrals  &
                              + memory_Sk_integral

      End Subroutine alloc_Rk_integrals


!====================================================================
      Subroutine alloc_Rk_integral(ns,ks)
!====================================================================
!     Allocates/deallocate space for B-spline Rk-integrals
!     for a single type and multipole index.
!--------------------------------------------------------------------
      Use DBS_integrals

      Implicit none
      Integer, intent(in) :: ns,ks

      if(associated(rkb)) Nullify(rkb)
      krk=-100; itype='aaaa'

      if(ntype1.ne.0) then
       Deallocate (rka1)
       krk1 = -100;  ntype1 = 0;  itype1='aaaa'
       memory_Rk_integral = 0.d0
      end if

      if(ns.gt.0) then
       Allocate(rka1(ns,ns,ks,ks))
       krk1 = -100;  ntype1 = 1;  itype1='bbbb'
       memory_Rk_integral = ns*ns*ks*ks*8/(1024d0*1024d0)
       Call alloc_DBS_moments
      end if

      memory_DBS_integrals =    memory_Rk_integrals  &
                              + memory_Rk_integral   &
                              + memory_Sk_integrals  &
                              + memory_Sk_integral

      End Subroutine alloc_Rk_integral


!====================================================================
      Subroutine alloc_Sk_integrals(ns,ks,kmin,kmax,ktype)
!====================================================================
!     Allocates/deallocate space for B-spline Rk-integrals
!     for different types and multipole indexes.
!     We expecting two types: 1(5). ppqq  and  2(6). pqqp
!--------------------------------------------------------------------
      Use DBS_integrals

      Implicit none
      Integer, intent(in) :: ns,ks,kmin,kmax,ktype

      if(associated(rkb)) Nullify(rkb)
      krk=-100
      itype='aaaa'

      if(allocated(isk)) then
       Deallocate (isk,rks)
       stype = 0
       memory_Sk_integrals = 0.d0
      end if

      if(ktype.gt.0) then
       Allocate(isk(kmin:kmax,ktype))
       Allocate(rks(ns,ns,2*ks-1,2*ks-1,kmin:kmax,ktype))
       krs_min = kmin
       krs_max = kmax
       isk = -100
       stype = ktype
       memory_Sk_integrals = ((kmax-kmin+1)*stype*4 + &
        ns*ns*(2*ks-1)*(2*ks-1)*(kmax-kmin+1)*stype*8)/(1024d0*1024d0)
       Call alloc_DBS_moments
      end if

      memory_DBS_integrals =    memory_Rk_integrals  &
                              + memory_Rk_integral   &
                              + memory_Sk_integrals  &
                              + memory_Sk_integral

      End Subroutine alloc_Sk_integrals



!====================================================================
      Subroutine alloc_Sk_integral(ns,ks)
!====================================================================
!     Allocates/deallocate space for B-spline Rk-integrals
!     for a single type and multipole index.
!--------------------------------------------------------------------
      Use DBS_integrals

      Implicit none
      Integer, intent(in) :: ns,ks

      if(associated(rkb)) Nullify(rkb)
      krk=-100; itype='aaaa'

      if(jtype1.ne.0) then
       Deallocate (rks1)
       krs1 = -100;  jtype1 = 0;  stype1='aaaa'
       memory_Sk_integral = 0.d0
      end if

      if(ns.gt.0) then
       Allocate(rks1(ns,ns,2*ks-1,2*ks-1))
       krs1 = -100;  jtype1 = 1;  stype1='bbbb'
       memory_Sk_integral = ns*ns*(2*ks-1)*(2*ks-1)*8/(1024d0*1024d0)
       Call alloc_DBS_moments
      end if

      memory_DBS_integrals =    memory_Rk_integrals  &
                              + memory_Rk_integral   &
                              + memory_Sk_integrals  &
                              + memory_Sk_integral

      End Subroutine alloc_Sk_integral


