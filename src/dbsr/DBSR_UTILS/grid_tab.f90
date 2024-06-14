!======================================================================
!     PROGRAM       G R I D _ T A B
!
!               C O P Y R I G H T -- 2024
!
!     Written by:   Stefanos Carlström
!                   email: stefanos.carlstrom@gmail.com
!
!======================================================================
!  Tabulate the basis functions of a grid
!======================================================================
!
!  INPUT FILES:
!
!     grid.dat     -  B-spline parameters of grid (knot.dat format)
!
!  OUTPUT FILES:
!
!     grid.dat.tab -  tabulated basis functions
!
!======================================================================
program grid_tab
  use DBS_grid
  use DBS_gauss
  use grid_tools

  Implicit none
  Integer, External :: Icheck_file

  character(256) :: grid_file="knot.dat", tab_file="",  name=""

  Real(8), allocatable :: r(:),cp(:),cq(:),chip(:,:),chiq(:,:)
  Real(8) :: tr, dr
  Integer :: nr = 101, i, j, nu, distribution=1, verbosity=0

  Real(8), external :: bvalu2

  call read_name(name)
  if(name == "?") call print_help

  Call Read_aarg('grid', grid_file)
  Call Read_iarg('nr', nr)
  Call Read_iarg('distribution', distribution)
  Call Read_iarg('verbosity', verbosity)

  if(grid_file == "") call print_help

  call check_file(grid_file)
  call load_grid(grid_file, verbosity)

  if(verbosity>1) then
     write(*,*) "grid_type = ", grid_type

     write(*,*) "ks = ", ks
     write(*,*) "ns = ", ns
     write(*,*) "ms = ", ms
     write(*,*) "nv = ", nv

     write(*,*) "ml = ", ml
     write(*,*) "me = ", me

     write(*,*) "hi = ", hi
     write(*,*) "he = ", he
     write(*,*) "hmax = ", hmax
     write(*,*) "rmax = ", rmax
     write(*,*) "tmax = ", tmax

     write(*,*) "ksp = ", ksp
     write(*,*) "nsp = ", nsp
     write(*,*) "ksq = ", ksq
     write(*,*) "nsq = ", nsq
  end if

  if(nr.le.0) stop "Number of steps must be larger than zero"

  allocate(r(nr),cp(nsp),cq(nsq),chip(nr,nsp),chiq(nr,nsq))

  call tabulation_grid(r, tmax, distribution)

  do j=1,nsp
     cp = 0
     cp(j) = 1
     do i=1,nr
        chip(i,j) = bvalu2(tp, cp, nsp, ksp, r(i), 0)
     end do
  end do

  do j=1,nsq
     cq = 0
     cq(j) = 1
     do i=1,nr
        chiq(i,j) = bvalu2(tq, cq, nsq, ksq, r(i), 0)
     end do
  end do

  nu = 1
  write(tab_file,'(a,".tab")') trim(grid_file)
  open(nu, file=trim(tab_file), form='FORMATTED')
  write(nu, '("nr = ",i5)') nr
  write(nu, '("nsp = ",i5)') nsp
  write(nu, '("nsq = ",i5)') nsq
  do i=1,nr
     write(nu,'(i5,1e26.13)') i, r(i)
     write(nu,'(4e26.13)') chip(i,:)
     write(nu,'(4e26.13)') chiq(i,:)
  end do
  close(nu)

  deallocate(r,cp,cq,chip,chiq)
contains
  subroutine print_help
    use grid_tools

    implicit none (external)

    write(*,*) "Usage: grid_tab grid= nr= distribution="
    write(*,*) ""
    call print_distribution_help
    stop
  end subroutine print_help

  subroutine load_grid(filename, verbosity)
    use DBS_grid
    use DBS_gauss

    implicit none
    Character(256), intent(in) :: filename
    Integer, intent(in) :: verbosity
    Character(256) :: name = ' '
    Real(8) :: z = 0.d0, awt = 0.d0

    if(verbosity>0) write(*,'("Loading knot set from ",a)') trim(filename)

    call def_grid(filename, name, z, awt)
    call alloc_DBS_gauss
    if(verbosity>2) then
       write(*,*) "Knot sequence:"
       write(*,'(4e26.13)') t
       write(*,*) "P knot sequence:"
       write(*,'(4e26.13)') tp
       write(*,*) "Q knot sequence:"
       write(*,'(4e26.13)') tq
    end if
  end subroutine load_grid
end program
