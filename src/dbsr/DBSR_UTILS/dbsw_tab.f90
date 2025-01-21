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
  use grid_tools

  Implicit none
  Character(1) :: ans
  Character(5) :: EL
  Character(256) :: AF="",BF="",knot_file="knot.dat"
  Real(8), allocatable ::  R(:),P(:),Q(:)
  Real(8) :: e_orb(1000)
  Real(8) :: S

  Real(8) :: tr, dr, rmax_tab
  Integer :: nr = 101, i, j, nu, distribution=1, verbosity=0

  Integer :: iout, nuw

  Real(8), external :: bvalu2

  ! ... input data:

  AF = '?'
  call Read_name(AF)
  write(*,*) "AF = ", AF

  if(AF.eq.'?') call print_help

  Call Read_aarg('grid', knot_file)

  ! ... set up B-splines:

  Call load_grid(knot_file)
  ! Call alloc_DBS_galerkin

  rmax_tab = tmax

  Call Read_iarg('nr', nr)
  Call Read_iarg('distribution', distribution)
  Call Read_rarg('rmax_tab', rmax_tab)
  Call Read_iarg('verbosity', verbosity)

  write(*,'("Tabulating orbitals for r = 0.0 .. ",f10.5)') rmax_tab

  ! ... radial w.f.:
  write(*,'("Loading orbitals from ",a)') trim(AF)
  nuw=1
  Open(nuw,file=trim(AF),status='OLD',form='UNFORMATTED')
  Call Read_pqbs(nuw, e_orb)
  Close(nuw)

  ! ... sets up grid points and initializes the values of the spline:
  allocate(R(NR),P(NR),Q(NR))
  call tabulation_grid(R, rmax_tab, distribution)

  ! ... Cycle over nl in input:

  BF=AF

  Do i=1,nbf
     do j=1,nr
        P(j) = bvalu2(tp, pq(1,1,i), nsp, ksp, R(j), 0)
        Q(j) = bvalu2(tq, pq(1,2,i), nsq, ksq, R(j), 0)
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
  subroutine print_help
    use grid_tools

    implicit none(external)

    write(*,*)
    write(*,*) 'dbsw_tab converts the B-spline radial orbital wbs-files into'
    write(*,*) 'text files suitable for graphic display'
    write(*,*)
    write(*,*) 'Call as:   dbsw_tab name.bsw  grid= nr= distribution= rmax_tab='
    write(*,*)
    write(*,*) 'file will be created for each orbital'
    write(*,*)
    write(*,*) 'rmax_tab is by default tmax, i.e. the end of the B-spline knot set'
    write(*,*) ""
    call print_distribution_help
    stop
  end subroutine print_help

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
