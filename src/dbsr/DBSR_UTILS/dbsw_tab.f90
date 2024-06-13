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
  Character(256) :: AF,BF,knot_file
  Real(8), allocatable ::  R(:),P(:),Q(:)
  Real(8) :: e_orb(1000)

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
  write(*,'("Loading orbitals from ",a)') trim(AF)
  nuw=1
  Open(nuw,file=trim(AF),status='OLD',form='UNFORMATTED')
  Call Read_pqbs(nuw, e_orb)
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
     write(*,*) trim(BF)
     iout=2; Open(iout,file=trim(BF))
     write(iout,'(3(5x,a,10x))') 'R','P','Q'

     S=1.d0
     !       S = c_au * 2.d0 * t(ns+1) / kbs(i)

     Do j=1,nr
        write(iout,'(3E26.16e3)') R(j),P(j),Q(j)*S
     End do

  End do


contains

  subroutine load_grid(filename)
    use DBS_grid
    use DBS_gauss

    implicit none
    Character(256), intent(in) :: filename
    Character(256) :: name = ' '
    Real(8) :: z = 0.d0, awt = 0.d0

    write(*,'("Loading knot set from ",a)') trim(filename)

    call Check_file(filename)
    call def_grid(filename, name, z, awt)
    call alloc_DBS_gauss

  end subroutine load_grid
end program dbsw_tab
