!======================================================================
!     UTILITY       D B S W _ T A B
!
!               C O P Y R I G H T -- 2008
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!     connverts the B-spline radial orbital wbs-files into tab-files
!     suitable for grafic display
!----------------------------------------------------------------------
!
!     INPUT FILE:    name.bsw
!
!     OUTPUT FILES:  name.bsw.nl  for each nl
!
!----------------------------------------------------------------------
!     ARGUMENTS:     name.bsw
!----------------------------------------------------------------------
program dbsw_tab
  Use DBS_grid
  Use DBS_gauss
  Use DBS_orbitals_pq
  Use zconst, only: c_au

  Implicit real(8) (A-H,O-Z)
  Character(1) :: ans
  Character(5) :: EL
  Character(40) :: AF,BF,knot_file
  Real(8), allocatable ::  R(:),P(:),Q(:)

  Integer i1,j1

  ! ... input data:

  AF = '?'
  call Read_name(AF)

  if(AF.eq.'?') then
     write(*,*)
     write(*,*) 'dbsw_tab connverts the B-spline radial orbital wbs-files into'
     write(*,*) 'text files suitable for grafic display'
     write(*,*)
     write(*,*) 'Call as:   dbsw_tab  name.bsw'
     write(*,*)
     write(*,*) 'file will be created for each orbital'
     Stop ' '
  end if

  ! ... set up B-splines:

  knot_file = ""
  Call Read_aarg('grid', knot_file)
  if(knot_file == "") then
     knot_file = 'knot.dat'
  end if

  Call load_grid(knot_file)
  ! Call alloc_DBS_galerkin

  ! ... radial w.f.:

  nuw=1
  Open(nuw,file=AF,status='OLD',form='UNFORMATTED')
  Call Read_pqbs(nuw)
  Close(nuw)

  ! ... sets up grid points and initializes the values of the spline:

  NR = nv*ks+2; Allocate(R(NR),P(NR),Q(NR))
  ii=1; R(1)=0.d0
  Do i=1,nv
     Do j=1,ks
        ii=ii+1
        R(ii) = gr(i,j)
     End do
  End do
  ii=ii+1; R(ii) = t(ns+1)

  ! ... Cycle over nl in input:

  BF=AF

  Do i=1,nbf
     ii = 1
     do i1=1,nv
        do j1=1,ks
           P(ii) = bvalu2(tp, pq(1,1,i), nsp, ksp, gr(i1,j1), 0)
           Q(ii) = bvalu2(tq, pq(1,2,i), nsq, ksq, gr(i1,j1), 0)
           ii = ii + 1
        end do
     end do


     write(BF,'(a,".",a)') trim(AF), adjustl(trim(ebs(i)))
     write(*,*) BF
     iout=2; Open(iout,file=BF)
     write(iout,'(3(5x,a,10x))') 'R','P','Q'

     S=1.d0
     !       S = c_au * 2.d0 * t(ns+1) / kbs(i)

     Do j=1,nr
        write(iout,'(3D16.8)') R(j),P(j),Q(j)*S
     End do

  End do


contains

  subroutine load_grid(filename)
    use DBS_grid
    use DBS_gauss

    implicit none
    Character(40), intent(in) :: filename
    Character(40) :: name = ' '
    Real(8) :: z = 0.d0, awt = 0.d0

    write(*,'("Loading knot set from ",a)') filename

    call Check_file(filename)
    call def_grid(filename, name, z, awt)
    call alloc_DBS_gauss

  end subroutine load_grid


  !======================================================================
  Subroutine Read_pqbs(nu)
    !======================================================================
    !
    !     read B-spline w.f. from bsw-file (unit nu) only those orbitals
    !     which labels are already in the list
    !
    !----------------------------------------------------------------------

    USE DBS_grid
    USE DBS_gauss
    USE DBS_orbitals_pq

    Implicit none

    Integer, intent(in) :: nu
    Integer :: i,j,k,l,n,m,itype,nsw,ksw,mw,kp,kq
    Character(5) :: elw
    Integer, External :: Ifind_bsorb,Iadd_bsorb
    Real(8) :: tt(ns+ks)

    rewind(nu)
    read(nu) itype,nsw,ksw,tt,kp,kq
    if(itype.ne.grid_type) Stop ' Read_pqbs:  another grid_type'
    if(ksw.ne.ks) Stop ' Read_pqbs:  ksw <> ks'
    if(nsw.ne.ns) Stop ' Read_pqbs:  nsw <> ns'
    if(ksp.ne.kp) Stop ' Read_pqbs:  ksp <> kp'
    if(ksq.ne.kq) Stop ' Read_pqbs:  ksq <> kq'
    Do i=1,ns+ks
       if(tt(i).ne.t(i)) Stop ' Read_pqbs:  t <> tt'
    End do

1   read(nu,end=2) elw,mw
    Call EL_NLJK(elw,n,k,l,j,i)
    m = Ifind_bsorb(n,k,i,2)
    mbs(m)=mw
    pq(1:ns,1,m)=0.d0; read(nu) pq(1:mw,1,m)
    pq(1:ns,2,m)=0.d0; read(nu) pq(1:mw,2,m)
    go to 1
2   Close(nu)
  End subroutine Read_pqbs
end program dbsw_tab
