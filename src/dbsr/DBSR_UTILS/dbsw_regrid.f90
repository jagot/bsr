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

  character(80) :: src_wfn_file="", dst_wfn_file="", &
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

  ! We first need to load the destination grid, to generate the
  ! quadrature nodes, on which we evaluate the orbitals.
  call load_grid(dst_grid_file)
  dst_nv = nv
  dst_ks = ks
  allocate(dst_gr(dst_nv,dst_ks))
  dst_gr(:,:) = gr(:,:)

  ! We then load the source grid, so we may evaluate the orbitals.
  call load_grid(src_grid_file)
  call load_source_orbitals(src_wfn_file, e_orb)
  norb = nbf
  call evaluate_orbitals(dst_gr, pr, qr, dst_nv, dst_ks)

  ! Finally, we load the destination grid again, to perform the
  ! inner products and solve the equation systems.
  call load_grid(dst_grid_file)
  call reallocate_orbitals(norb, ns)
  call project_orbitals(pr, qr)
  deallocate(dst_gr,pr,qr)

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
    Character(40), intent(in) :: filename
    Character(40) :: name = ' '
    Real(8) :: z = 0.d0, awt = 0.d0

    write(*,'("Loading knot set from ",a)') filename

    call def_grid(filename, name, z, awt)
    call alloc_DBS_gauss

    ! write(*,*) "Knot set:"
    ! write(*,'(5(1F26.14))') t
    ! ! write(*,*) "Quadrature points:"
    ! ! write(*,'(5(1F26.14))') gr
    ! ! write(*,*) "Quadrature weights:"
    ! ! write(*,'(5(1F26.14))') grw

  end subroutine load_grid

  subroutine load_source_orbitals(filename, e_orb)
    use DBS_orbitals_pq

    Implicit none
    character(80), intent(in) :: filename
    real(8), intent(out) :: e_orb(:)

    Integer :: nu = 1

    open(nu, file=trim(filename), form='UNFORMATTED', status='OLD')
    call read_pqbs(nu, e_orb)
    close(nu)

    write(*,'("Number of orbitals: ",i10)') nbf
  end subroutine load_source_orbitals

  !======================================================================
  Subroutine read_pqbs(nu, e_orb)
    !======================================================================
    !     read B-spline w.f. from bsw-file (unit nu)
    !     stolen from bsw_rw.f90
    !----------------------------------------------------------------------
    Use DBS_grid
    Use DBS_gauss
    Use DBS_orbitals_pq

    Implicit none
    Integer, intent(in) :: nu
    real(8), intent(out) :: e_orb(:)

    Integer :: i,j,k,l,n,m,itype,nsw,ksw,mw,kp,kq
    Character(5) :: elw
    Integer, external :: Ifind_bsorb
    Real(8) :: tt(ns+ks), S

    rewind(nu)
    read(nu) itype,nsw,ksw,tt,kp,kq
    if(grid_type.gt.0.and.itype.ne.grid_type) &
         Stop 'Stop in read_dbsw: another grid_type ?'
    if(ksw.ne.ks) Stop ' Read_pqbs:  ksw <> ks'
    if(nsw.ne.ns) Stop ' Read_pqbs:  nsw <> ns'
    if(ksp.ne.kp) Stop ' Read_pqbs:  ksp <> kp'
    if(ksq.ne.kq) Stop ' Read_pqbs:  ksq <> kq'
    k=1
    Do i=1,ns+ks
       if(abs(t(i)-tt(i)).lt.1.d-12) Cycle; k=0; Exit
    End do
    if(k.eq.0) Stop 'Stop in read_pqbs: another knot grid ?'

1   read(nu,end=2) elw,mw,S
    Call EL_NLJK(elw,n,k,l,j,i)
    m = Ifind_bsorb(n,k,i,2)
    e_orb(m) = S
    mbs(m)=mw
    pq(1:ns,1,m)=0.d0; read(nu) pq(1:mw,1,m)
    pq(1:ns,2,m)=0.d0; read(nu) pq(1:mw,2,m)
    bpq(:,1,m) = MATMUL(fpbs,pq(:,1,m))
    bpq(:,2,m) = MATMUL(fqbs,pq(:,2,m))
    go to 1
2   Close(nu)

  End subroutine read_pqbs

  subroutine evaluate_orbitals(dst_gr, pr, qr, dst_nv, dst_ks)
    Use DBS_grid
    Use DBS_gauss
    Use DBS_orbitals_pq

    Implicit none

    Real(8), allocatable, intent(in) :: dst_gr(:,:)
    Real(8), allocatable, intent(out) :: pr(:,:,:), qr(:,:,:)
    Integer, intent(in) :: dst_nv, dst_ks

    Integer :: io, i, j, m

    Real(8), external :: bvalu2

    allocate(pr(dst_nv,dst_ks,nbf),qr(dst_nv,dst_ks,nbf))

    pr = 0.d0
    qr = 0.d0

    write(*,*) "nbf = ", nbf, "mbf = ", mbf
    do io=1,nbf
       m = mbs(io)
       do i=1,dst_nv
          do j=1,dst_ks
             pr(i,j,io) = bvalu2(tp, pq(1,1,io), nsp, ksp, dst_gr(i,j), 0)
             qr(i,j,io) = bvalu2(tq, pq(1,2,io), nsq, ksq, dst_gr(i,j), 0)
          end do
       end do
    end do
  end subroutine evaluate_orbitals

  subroutine reallocate_orbitals(norb, ns)
    Use DBS_orbitals_pq

    Implicit none

    Integer, intent(in) :: norb, ns

    Integer, allocatable :: tmp_mbs(:)      ! number of splines
    Character(5), allocatable :: tmp_ebs(:) ! spectroscopic notation

    allocate(tmp_mbs(norb), tmp_ebs(norb))
    tmp_mbs(:) = mbs(1:norb)
    tmp_ebs(:) = ebs(1:norb)

    call alloc_DBS_orbitals_pq(0, ns)
    call alloc_DBS_orbitals_pq(norb, ns)
    nbf = norb
    mbs(1:norb) = ns
    ebs(1:norb) = tmp_ebs(:)

    deallocate(tmp_mbs, tmp_ebs)
  end subroutine reallocate_orbitals

  subroutine project_orbitals(pr, qr)
    Use DBS_grid
    Use DBS_gauss
    Use DBS_orbitals_pq

    Implicit none

    Real(8), allocatable, intent(in) :: pr(:,:,:), qr(:,:,:)

    Integer :: io, i, j, iv, ip

    Real(8), external :: bvalu2, QUADR

    Real(8) :: Pcoef(ns), Qcoef(ns), a(ns,ns)
    Real(8) :: pm, qm, d, PP, QQ, PN

    do io = 1,nbf
       Pcoef = 0.d0
       Qcoef = 0.d0
       do iv = 1,nv
          do ip = 1,ksp
             i = iv+ip-1
             do j = 1,ks
                Pcoef(i) = Pcoef(i) + pr(iv,j,io)*grw(iv,j)*pbsp(iv,j,ip)
                Qcoef(i) = Qcoef(i) + qr(iv,j,io)*grw(iv,j)*qbsp(iv,j,ip)
             end do
          end do
       end do

       a(1:nsp-1,1:nsp-1)=fpbs(2:nsp,2:nsp)
       Call gaussj (a,nsp-1,ns,Pcoef(2),1,ns)
       Pcoef(1)=0.d0

       a(1:nsq-1,1:nsq-1)=fqbs(2:nsq,2:nsq)
       Call gaussj (a,nsq-1,ns,Qcoef(2),1,ns)
       Qcoef(1)=0.d0

       pq(:,1,io) = Pcoef
       pq(:,2,io) = Qcoef

       ! ... check the normalization

       PN = sqrt(QUADR(pq(1,1,io),pq(1,1,io),0))
       pq(:,:,io)=pq(:,:,io)/PN

       ! ... max.deviation from original orbital:

       pm = 0.d0
       do i = 1,nv
          do j = 1,ks
             PP = bvalu2(tp,pq(1,1,io),nsp,ksp,gr(i,j),0)
             d = abs(PP-pr(i,j,io));  if(d.gt.pm) pm = d

             QQ = bvalu2(tq,pq(1,2,io),nsq,ksq,gr(i,j),0)
             d = abs(QQ-qr(i,j,io)); if(d.gt.qm) qm = d
          end do
       end do

       write(*,'(a5,4(a,E12.3))')  &
            ebs(io), '  diff_p =',pm,'   diff_q =',qm,'   (norm -1) =',PN-1.d0
    end do
  end subroutine project_orbitals

  subroutine save_orbitals(filename, e_orb)
    Use DBS_grid
    Use DBS_orbitals_pq

    Implicit none

    Character(40), intent(in) :: filename
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
