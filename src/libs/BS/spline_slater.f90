!====================================================================
    MODULE spline_slater
!====================================================================
!   contains arrays and variables needed for calculations of 
!   Breit-Pauli integrals by direct integration in B-spline basis
!--------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: ic1=0,ic2=0, iy1=0,iy2=0, ky=-100, kz=-100, ktx = 0

! . the above integers serve to check the repeated calculations for
! . arrays fc1,fc2 (ic1,ic2); fy1,fy2 (iy1,iy2), dyk,ipvtd (ky),
! . dz,ipvtz (kz);  ktx = 3*ks-2.

    REAL(8), ALLOCATABLE :: fc1(:,:),fc2(:,:), fy1(:,:),fy2(:,:)

!   fc1,fc2, fy1,fy2 - serve for storing the gauss-points (nv,ks)
!                      representation of radial orbitals

    REAL(8), ALLOCATABLE ::  fc(:,:)      ! work array

    REAL(8), ALLOCATABLE :: fyk(:,:), yk(:)

!   fyk - the gauss-point represantation of the Yk or Zk functions
!   yk  - the B-spline represantation of the Yk or Zk functions

    REAL(8), ALLOCATABLE :: dyk(:,:), dzk(:,:)
    INTEGER, ALLOCATABLE :: ipvtd(:), ipvtz(:)

!   dyk(ktx,ns),ipvtd(ns) - factorization of the Yk operator:
!                           Yk --> d^2/dr^2 - k(k+1)/r^2

!   dzk(ktx,ns),ipvtz(ns) - factorization of the Zk operator:
!                           Yk --> d/dr + k/r

    END MODULE spline_slater


!====================================================================
      SUBROUTINE allocate_slater
!====================================================================
!     allocates (reallocate) space for arrais in MODULE spline_slater
!--------------------------------------------------------------------
      USE spline_param
      USE spline_slater

      ktx = 3*ks-2
      if(Allocated(fc1)) Deallocate(fc1,fc2,fy1,fy2,fc, yk,fyk, &
                                    dyk,ipvtd, dzk,ipvtz)
      ALLOCATE(fc1(nv,ks), fc2(nv,ks), fy1(nv,ks), fy2(nv,ks),  &
               fc(nv,ks), yk(ns), fyk(nv,ks),                   &
               dyk(ktx,ns),ipvtd(ns),dzk(ktx,ns),ipvtz(ns))
	  fc1 = 0.d0; fc2 = 0.d0; fy1 = 0.d0; fy2 = 0.d0
	  fc = 0.d0; yk = 0.d0; fyk = 0.d0		   
      ic1=0; ic2=0; iy1=0; iy2=0; ky=-100; kz=-100
      
      END SUBROUTINE allocate_slater


!====================================================================
      SUBROUTINE dealloc_slater
!====================================================================
!     deallocates space for arrais in MODULE spline_slater
!--------------------------------------------------------------------
      USE spline_slater

      if(allocated(fc1)) Deallocate(fc1,fc2,fy1,fy2,fc, yk,fyk, &
                                    dyk,ipvtd, dzk,ipvtz)

      ic1=0; ic2=0; iy1=0; iy2=0; ky=-100; kz=-100; ktx = 0
      
      END SUBROUTINE dealloc_slater

