!======================================================================
!     PROGRAM       D B S W _ R E G R I D
!
!               C O P Y R I G H T -- 2024
!
!     Written by:   Stefanos Carlström
!                   email: stefanos.carlstrom@gmail.com
!
!======================================================================
!  Transform a wavefunction onto a different grid
!======================================================================
!
!  INPUT FILES:
!
!     src_wfn  -  w.f. file on source grid
!     src_grid -  B-spline parameters of source grid
!     dst_grid -  B-spline parameters of destination grid
!
!  OUTPUT FILES:
!
!     new.bsw     -  w.f. on destination grid
!
!======================================================================
program dbsw_grid
  use DBS_grid
  use DBS_gauss
  use DBS_orbitals_pq

  Implicit none
  Integer, External :: Icheck_file

  character(256) :: src_wfn_file="", dst_wfn_file="", &
       src_grid_file="", dst_grid_file="", &
       name=""

  Real(8), allocatable :: dst_gr(:,:), pr(:,:,:), qr(:,:,:)
  Real(8) :: e_orb(1000)
  Integer :: dst_nv, dst_ks, norb

  call read_name(name)
  if(name == "?") call print_help

  Call Read_aarg('src_wfn', src_wfn_file)
  Call Read_aarg('dst_wfn', dst_wfn_file)
  Call Read_aarg('src_grid', src_grid_file)
  Call Read_aarg('dst_grid', dst_grid_file)

  if(src_wfn_file == "" .or. dst_wfn_file == "" .or. &
       src_grid_file == "" .or. dst_grid_file == "") then
     call print_help
  end if

  call check_file(src_wfn_file)
  call check_file(src_grid_file)
  call check_file(dst_grid_file)

  ! if(Icheck_file(dst_wfn_file) == 1) then
  !    write(*,'("Destination wfn already exists at ",a)') dst_wfn_file
  !    write(*,*) "Please delete the file, if you wish to regrid again"
  !    stop 1
  ! end if

  call load_grid(dst_grid_file)
  call load_source_orbitals(src_wfn_file, e_orb)

  ! Save the regridded orbitals to file
  call save_orbitals(dst_wfn_file, e_orb)

contains

  subroutine print_help
    write(*,*) "Usage: dbsw_regrid src_wfn= dst_wfn= src_grid= dst_grid="
    stop ''
  end subroutine print_help

  subroutine load_grid(filename)
    use DBS_grid
    use DBS_gauss

    implicit none
    Character(256), intent(in) :: filename
    Character(40) :: name = ' '
    Real(8) :: z = 0.d0, awt = 0.d0

    write(*,'("Loading knot set from ",a)') trim(filename)

    call def_grid(filename, name, z, awt)
    call alloc_DBS_gauss
  end subroutine load_grid

  subroutine load_source_orbitals(filename, e_orb)
    use DBS_orbitals_pq

    Implicit none
    character(256), intent(in) :: filename
    real(8), intent(out) :: e_orb(:)

    Integer :: nu = 1

    open(nu, file=trim(filename), form='UNFORMATTED', status='OLD')
    call read_pqbs(nu, e_orb, .true.)
    close(nu)

    write(*,'("Number of orbitals: ",i10)') nbf
  end subroutine load_source_orbitals

  subroutine save_orbitals(filename, e_orb)
    Use DBS_grid
    Use DBS_orbitals_pq

    Implicit none

    Character(256), intent(in) :: filename
    real(8), intent(in) :: e_orb(:)

    Integer :: nu=1, i

    open(nu, file=trim(filename), form='UNFORMATTED')
    write(nu) grid_type,ns,ks,t(1:ns+ks),ksp,ksq
    do i=1,nbf
       write(nu) ebs(i),mbs(i),e_orb(i)
       write(nu) pq(1:mbs(i),1,i)
       write(nu) pq(1:mbs(i),2,i)
    end do
    close(nu)
  end subroutine save_orbitals
end program dbsw_grid
