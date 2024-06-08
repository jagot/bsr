!======================================================================
!     UTILITY    B S W _ R W
!
!               C O P Y R I G H T -- 2020
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!======================================================================
!     Rewrite the  double (p,q) B-spline representation for orbitals
!     in GRASP package format, w-files
!----------------------------------------------------------------------
!
!     INPUT ARGUMENTS:
!
!     name.bsw              -  input w-file (or file with list of w-files)
!
!     INPUT FILES:
!
!     name.bsw              -  B-spline representation for orbitals
!     knot.dat              -  parameters of B-splines
!
!     OUTPUT FILES:
!
!     name.w                -  w-files
!
!---------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_orbitals_pq

      Implicit real(8) (A-H,O-Z)

      Real(8) :: e(1000) ! ???

! ... files:

      Character(80) :: A
      Integer :: nub =1;  Character(80) :: AF_bsw = 'name.bsw'
      Integer :: nuw =2;  Character(80) :: AF_w   = 'name.w'

!----------------------------------------------------------------------
! ... check input data:

      Call Read_name(A)

      if(A.eq.'?') then
        write(*,*)
        write(*,*) 'bsw_rw: convert the DBSR bsw-file to GRASP w-format'
        write(*,*)
        write(*,*) 'knot.dat is needed'
        write(*,*)
        write(*,*) 'Call as:  bsw_rw  name.bsw'
        write(*,*)
        write(*,*) 'OUTPUT:  name.w '
        Stop ' '
       end if

!----------------------------------------------------------------------
! ... sets up grid points and initializes the values of the spline:

      Call Read_knot_dat

      Call alloc_DBS_gauss

!----------------------------------------------------------------------
! ... read arguments and define fails:


      AF_bsw=trim(A);  i = LEN_TRIM(A)
      AF_w = A(1:i-4)//'.w'

      Open(nub,file=AF_bsw,form='UNFORMATTED',status='OLD')

! ... read w.f. in B-spline format:

      Call Read_pqbs(nub,e)

      Open(nuw,file=AF_w,form='UNFORMATTED')

! ... output in the GRASP format:

      Call GRASP_wfn(nuw,e)

      End  !  program  pqbs_rw

!======================================================================
      Subroutine GRASP_wfn(nuw,e)
!======================================================================
! ... output in GRASP format:
!----------------------------------------------------------------------
      Use zconst,      only: c_au
      Use DBS_nuclear, only: nuclear, atomic_number
      Use DBS_grid
      Use DBS_gauss
      Use DBS_orbitals_pq

      Implicit none

      Integer, parameter :: ng = 540  ! max. number of points in GRASP
      Real(8) :: yp(ng),yq(ng),r(ng), e(*)
      Real(8) :: P0, gamma, r_max, RNT,HNT, z
      Integer :: i,j, io,m,np,nr, nuw
      Real(8), external :: bvalu2

! ... radial points for output:

      z = atomic_number
      if(nuclear.eq.'point') then
       RNT = EXP (-65.0d0/16.0d0) / z
       HNT = 0.5d0**4
       np  = ng
      else
       RNT = 2.d-6
       HNT = 5.d-2
       np  = min(ng,220)
      end if

      r_max = 0.d0
      Do io = 1,nbf
       if(r_max.lt.t(mbs(io)+ks)) r_max=t(mbs(io)+ks)
      End do
      Do
       Do i=1,np
        r(i)=RNT*(exp((i-1)*HNT)-1.d0)
       End do
       if(r(np).gt.r_max) Exit
       HNT = HNT*1.2
      End do
      r(1) = 0.d0

! ... output to unit nuw:

      write(nuw) 'G92RWF'

! ... cycle over orbitals

      Do io = 1,nbf; m = mbs(io)

       r_max = t(m+ks)
write(*,*) ebs(io), m, r_max

       yp = 0.d0; if(m.lt.ns) pq(m+1:ns,1,io)=0.d0
       yq = 0.d0; if(m.lt.ns) pq(m+1:ns,2,io)=0.d0
       Do i = 2,np-1
        yp(i) =  bvalu2 (tp, pq(1,1,io), nsp, ksp, r(i), 0)
        yq(i) =  bvalu2 (tq, pq(1,2,io), nsq, ksq, r(i), 0)
        nr = i; if(r(i+1).gt.r_max) Exit
       End do
       nr = nr + 1

       gamma = lbs(io) + 1
       if(nuclear.eq.'point') gamma = sqrt (kbs(io)**2 - (z/c_au)**2)
       P0 = yp(2)/r(2)**gamma

       write(nuw) nbs(io),kbs(io),e(io),nr
       write(nuw) P0,yp(1:nr),yq(1:nr)
       write(nuw) r(1:nr)

      End do


      End Subroutine GRASP_wfn
