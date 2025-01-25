module DBS_convolutions
  Implicit none
  Public convol
contains
  !======================================================================
  Subroutine convol(ns,ks,a,d,icase,sym_i,sym_j)
    !======================================================================
    !     convolutes the rkb(i,j,i',j') array of spline integrals
    !     with density matrix d(:,:)
    !
    !     results in array a(:,:)
    !
    !     icase =  1  - convolution other 2 and 4 variables, RK(.a;.b)
    !              2  - convolution other 1 and 3 variables, RK(a.;b.)
    !              3  - convolution other 2 and 3 variables, RK(.a;b.)
    !              4  - convolution other 1 and 4 variables, RK(a.;.b)
    !
    !     sym_i  ->  symmetry in respect of i,i' variables ('s','l','n')
    !     sym_j  ->  symmetry in respect of j,j' variables ('s','l','n')
    !
    !     combination of sym_i and sym_j leads to different represantation
    !     for a and d:   a(ns,ks),  a(ns,2*ks+1),  a(ns,ns)
    !----------------------------------------------------------------------
    Use DBS_integrals, only: rkb, itype, d24_rk
    Use DBS_debug

#ifdef DEBUG_SPEEDUPS
    Use Timer
#endif

    Implicit none
    Integer, intent(in) :: ns,ks,icase
    Character, intent(in) :: sym_i,sym_j
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)
    ! local variables
    Integer :: i,j, ip,jp, imin,imax, jmin,jmax, ii,jj
    Real(8) :: c,t1,t2
    Integer :: num_blas_threads

#ifdef DEBUG_SPEEDUPS
    Real(8), allocatable :: a_ref(:,:)
    Real(8) :: discrepancy
    Real(8), parameter :: tolerance = sqrt(epsilon(1.d0))
    Character(256) :: timer_label

    write(timer_label, '("convol c=",i0," si=",a," sj=",a)') &
         icase, sym_i, sym_j
    call TimerStart(trim(timer_label))
#endif

    Call CPU_time(t1)

    if(icase.le.2) a(1:ns,1:ks)=0.d0
    if(icase.gt.2) a(1:ns,1:ns)=0.d0

    Select case(icase)
       !----------------------------------------------------------------------
    Case(1)                                          !  I( . a ; . b)
       ! write(*,'("convol: I(.a;.b) ",a," ",a,", ks = ",i0,", ns = ",i0,", itype = ",a)') &
       !      sym_i, sym_j, ks, ns, itype

       if(sym_i.eq.'s'.and.sym_j.eq.'s') then
          call Ipapb_ss(ns, ks, a, d)
       elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then
          call Ipapb_nn(ns, ks, a, d)
       elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then
          call Ipapb_ll(ns, ks, a, d)
       elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then
          call Ipapb_ln(ns, ks, a, d)
       elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then
          call Ipapb_nl(ns, ks, a, d)
       end if
       !----------------------------------------------------------------------
    Case(2)                                         !  I( a . ; b . )
       ! write(*,'("convol: I(a.;b.) ",a," ",a,", ks = ",i0,", ns = ",i0,", itype = ",a)') &
       !      sym_i, sym_j, ks, ns, itype

       if(sym_i.eq.'s'.and.sym_j.eq.'s') then
          call Iapbp_ss(ns, ks, a, d)
       elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then
          call Iapbp_nn(ns, ks, a, d)
       elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then
          call Iapbp_ll(ns, ks, a, d)
       elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then
          call Iapbp_ln(ns, ks, a, d)
       elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then
          call Iapbp_nl(ns, ks, a, d)
       end if

       !----------------------------------------------------------------------
    Case(3)                                         ! I( . a ; b . )
       ! write(*,'("convol: I(.a;b.) ",a," ",a,", ks = ",i0,", ns = ",i0,", itype = ",a)') &
       !      sym_i, sym_j, ks, ns, itype
       a(1:ns,1:ns) = 0.d0
       if(sym_i.eq.'s'.and.sym_j.eq.'s') then
          call Ipabp_ss(ns, ks, a, d)
       elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then
          call Ipabp_ll(ns, ks, a, d)
       elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then
          call Ipabp_nn(ns, ks, a, d)
       elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then
          call Ipabp_nl(ns, ks, a, d)
       elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then
          call Ipabp_ln(ns, ks, a, d)
       end if
       !----------------------------------------------------------------------
    Case(4)                                         ! I( a . ; . b )
       ! write(*,'("convol: I(a.;.b) ",a," ",a,", ks = ",i0,", ns = ",i0,", itype = ",a)') &
       !      sym_i, sym_j, ks, ns, itype
       a(1:ns,1:ns) = 0.d0

       if(sym_i.eq.'s'.and.sym_j.eq.'s') then
          call Iappb_ss(ns, ks, a, d)
       elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then
          call Iappb_nn(ns, ks, a, d)
       elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then
          call Iappb_ll(ns, ks, a, d)
       elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then
          call Iappb_nl(ns, ks, a, d)
       elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then
          call Iappb_ln(ns, ks, a, d)
       end if
       !----------------------------------------------------------------------
    Case default

       Error Stop 'convol: unknown case'

    End Select

    Call CPU_time(t2)
    ic_convol = ic_convol + 1
    if(icase.le.2) time_convol=time_convol + (t2-t1)

#ifdef DEBUG_SPEEDUPS
    call TimerStop(trim(timer_label))
#endif
  End Subroutine convol

  ! ----------------------------------------------------------------------------------------------------
  ! Individual convolutions follow
  ! ----------------------------------------------------------------------------------------------------
  ! I(. a ; . b)
  ! ----------------------------------------------------------------------------------------------------
  subroutine Ipapb_ss(ns, ks, a, d)
    Use DBS_integrals, only: rkb, d24_rk
    Use BLAS
    ! Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)

    Integer :: i,ip,jp
    Integer :: num_blas_threads

    ! call TimerStart('convol: I(.a;.b) s s, step 1')

    a(1:ns,1:ks)=0.d0

    ! Compared to the old implementation, we do not limit the
    ! contraction over i to 1,ns-ip+1, where ip = 1,ks then has
    ! to be the outer loop. Instead we use the fact that the
    ! array storing the Rk-integrals, rkb, is zero outside that
    ! boundary, and let i run freely all the way to ns.

    num_blas_threads = blas_get_num_threads()
    call blas_set_num_threads(1)
    !$omp parallel do collapse(2) default(none) private(ip,jp) shared(ks,ns,rkb,d24_rk)
    do ip=1,ks
       do jp=1,ks
          d24_rk(:,ip,jp) = matmul(rkb(1:ns,1:ns,ip,jp), d(1:ns,jp))
       end do
    end do
    !$omp end parallel do
    call blas_set_num_threads(num_blas_threads)

    ! call TimerStop('convol: I(.a;.b) s s, step 1')
    ! call TimerStart('convol: I(.a;.b) s s, step 2')

    !$omp parallel do default(none) private(ip,i) shared(ks,ns,a,d24_rk)
    do ip=1,ks
       do i=1,ns-ip+1
          a(i,ip) = sum(d24_rk(i,ip,:))
       end do
    end do
    !$omp end parallel do

    ! call TimerStop('convol: I(.a;.b) s s, step 2')

#ifdef DEBUG_SPEEDUPS_CONVOL_PAPB_SS
    call Ipapb_ss_old(ns, ks, a, d)
#endif

  end subroutine Ipapb_ss

#ifdef DEBUG_SPEEDUPS_CONVOL_PAPB_SS
  subroutine Ipapb_ss_old(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(in) :: a(ns,*)

    Integer :: i,ip,jp

    Real(8), allocatable :: a_ref(:,:)
    Real(8) :: discrepancy
    Real(8), parameter :: tolerance = sqrt(epsilon(1.d0))

    call TimerStart('convol: I(.a;.b) s s, old')
    allocate(a_ref(1:ns,1:ks))
    a_ref = 0.d0
    do ip=1,ks
       do i=1,ns-ip+1
          a_ref(i,ip) = SUM(d(1:ns,1:ks)*rkb(i,1:ns,ip,1:ks))
       end do
    end do
    call TimerStop('convol: I(.a;.b) s s, old')
    discrepancy = sum(abs(a(1:ns,1:ks) - a_ref(1:ns,1:ks)))
    write(*,'("Discrepancy: ",e26.16,", tolerance: ",e26.16)') discrepancy, tolerance
    if(discrepancy > tolerance) error stop "Discrepancy exceeds tolerance"
    deallocate(a_ref)
  end subroutine Ipapb_ss_old
#endif
  ! --------------------------------------------------------------------------------
  subroutine Ipapb_nn(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)

    Integer :: i,j, ip,jp, imin,imax, jmin,jmax
    Real(8) :: c

    call TimerStart('convol: I(.a;.b) n n')
    a(1:ns,1:ks)=0.d0
    !$omp parallel do default(none) private(ip,imin,imax,i,c,jp,jmin,jmax) shared(ns,ks,a,d,rkb)
    do ip=1,ks+ks-1
       imin=max(1,1+ks-ip)
       imax=min(ns,ns+ks-ip)
       do i =imin,imax
          c=0.d0
          do jp=1,ks+ks-1
             jmin=max(1,1+ks-jp)
             jmax=min(ns,ns+ks-jp)
             c = c + sum(d(jmin:jmax,jp)*rkb(i,jmin:jmax,ip,jp))
          end do
          a(i,ip)=c
       end do
    end do
    !$omp end parallel do
    call TimerStop('convol: I(.a;.b) n n')

#ifdef DEBUG_SPEEDUPS
    call Ipapb_nn_old(ns, ks, a, d)
#endif
  end subroutine Ipapb_nn

#ifdef DEBUG_SPEEDUPS
  subroutine Ipapb_nn_old(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(in) :: a(ns,*)

    Integer :: i,j, ip,jp, imin,imax, jmin,jmax
    Real(8) :: c

    Real(8), allocatable :: a_ref(:,:)
    Real(8) :: discrepancy
    Real(8), parameter :: tolerance = sqrt(epsilon(1.d0))

    call TimerStart('convol: I(.a;.b) n n, old')
    allocate(a_ref(1:ns,1:ks))
    a_ref(1:ns,1:ks) = 0.d0
    do ip=1,ks+ks-1
       imin=max(1,1+ks-ip)
       imax=min(ns,ns+ks-ip)
       do i =imin,imax
          c=0.d0
          do jp=1,ks+ks-1
             jmin=max(1,1+ks-jp)
             jmax=min(ns,ns+ks-jp)
             do j =jmin,jmax
                c=c+d(j,jp)*rkb(i,j,ip,jp)
             end do
          end do
          a_ref(i,ip)=c
       end do
    end do
    call TimerStop('convol: I(.a;.b) n n, old')
    discrepancy = sum(abs(a(1:ns,1:ks) - a_ref(1:ns,1:ks)))
    write(*,'("Discrepancy: ",e26.16,", tolerance: ",e26.16)') discrepancy, tolerance
    if(discrepancy > tolerance) error stop "Discrepancy exceeds tolerance"
    deallocate(a_ref)
  end subroutine Ipapb_nn_old
#endif
  ! --------------------------------------------------------------------------------
  subroutine Ipapb_ll(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)

    Integer :: i,j, ip,jp
    Real(8) :: c

    call TimerStart('convol: I(.a;.b) l l')
    a(1:ns,1:ks)=0.d0
    do ip=1,ks
       do i=ks+1-ip,ns
          c=0.d0
          do jp=1,ks
             do j=ks+1-jp,ns
                c=c+d(j,jp)*rkb(i,j,ip,jp)
             end do
          end do
          a(i,ip)=c
       end do
    end do
    call TimerStop('convol: I(.a;.b) l l')
  end subroutine Ipapb_ll
  ! --------------------------------------------------------------------------------
  subroutine Ipapb_ln(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)

    Integer :: i,j, ip,jp, jmin,jmax
    Real(8) :: c

    call TimerStart('convol: I(.a;.b) l n')
    a(1:ns,1:ks)=0.d0
    do ip=1,ks
       do i=ks+1-ip,ns
          c=0.d0
          do jp=1,ks+ks-1
             jmin=max(1,1+ks-jp)
             jmax=min(ns,ns+ks-jp)
             do j =jmin,jmax
                c=c+d(j,jp)*rkb(i,j,ip,jp)
             end do
          end do
          a(i,ip)=c
       end do
    end do
    call TimerStop('convol: I(.a;.b) l n')
  end subroutine Ipapb_ln
  ! --------------------------------------------------------------------------------
  subroutine Ipapb_nl(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)

    Integer :: i,j, ip,jp, imin,imax
    Real(8) :: c

    call TimerStart('convol: I(.a;.b) n l')
    a(1:ns,1:ks)=0.d0
    do ip=1,ks+ks-1
       imin=max(1,1+ks-ip)
       imax=min(ns,ns+ks-ip)
       do i =imin,imax
          c=0.d0
          do jp=1,ks
             do j=ks+1-jp,ns
                c=c+d(j,jp)*rkb(i,j,ip,jp)
             end do
          end do
          a(i,ip)=c
       end do
    end do
    call TimerStop('convol: I(.a;.b) n l')
  end subroutine Ipapb_nl

  ! ----------------------------------------------------------------------------------------------------
  ! I(a . ; b .)
  ! ----------------------------------------------------------------------------------------------------
  subroutine Iapbp_ss(ns, ks, a, d)
    Use DBS_integrals, only: rkb, d24_rk
    Use BLAS
    ! Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)

    Integer :: i,ip,jp
    Integer :: num_blas_threads

    ! call TimerStart('convol: I(a.;b.) s s, step 1')

    a(1:ns,1:ks)=0.d0

    ! Compared to the old implementation, we do not limit the
    ! contraction over i to 1,ns-ip+1, where ip = 1,ks then has
    ! to be the outer loop. Instead we use the fact that the
    ! array storing the Rk-integrals, rkb, is zero outside that
    ! boundary, and let i run freely all the way to ns.

    num_blas_threads = blas_get_num_threads()
    call blas_set_num_threads(1)
    !$omp parallel do collapse(2) default(none) private(ip,jp) shared(ks,ns,rkb,d24_rk)
    do ip=1,ks
       do jp=1,ks
          d24_rk(:,ip,jp) = matmul(transpose(rkb(1:ns,1:ns,ip,jp)), d(1:ns,ip))
       end do
    end do
    !$omp end parallel do
    call blas_set_num_threads(num_blas_threads)

    ! call TimerStop('convol: I(a.;b.) s s, step 1')
    ! call TimerStart('convol: I(a.;b.) s s, step 2')

    !$omp parallel do default(none) private(ip,i) shared(ks,ns,a,d24_rk)
    do ip=1,ks
       do i=1,ns-ip+1
          a(i,ip) = sum(d24_rk(i,:,ip))
       end do
    end do
    !$omp end parallel do

    ! call TimerStop('convol: I(a.;b.) s s, step 2')

#ifdef DEBUG_SPEEDUPS_CONVOL_APBP_SS
    call Iapbp_ss_old(ns, ks, a, d)
#endif

  end subroutine Iapbp_ss

#ifdef DEBUG_SPEEDUPS_CONVOL_APBP_SS
  subroutine Iapbp_ss_old(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(in) :: a(ns,*)

    Integer :: i,ip,jp

    Real(8), allocatable :: a_ref(:,:)
    Real(8) :: discrepancy
    Real(8), parameter :: tolerance = sqrt(epsilon(1.d0))

    call TimerStart('convol: I(a.;b.) s s, old')
    allocate(a_ref(1:ns,1:ks))
    a_ref = 0.d0
    do ip=1,ks
       do i=1,ns-ip+1
          a_ref(i,ip) = SUM(d(1:ns,1:ks)*rkb(1:ns,i,1:ks,ip))
       end do
    end do
    call TimerStop('convol: I(a.;b.) s s, old')
    discrepancy = sum(abs(a(1:ns,1:ks) - a_ref(1:ns,1:ks)))
    write(*,'("Discrepancy: ",e26.16,", tolerance: ",e26.16)') discrepancy, tolerance
    if(discrepancy > tolerance) error stop "Discrepancy exceeds tolerance"
    deallocate(a_ref)
  end subroutine Iapbp_ss_old
#endif
  ! --------------------------------------------------------------------------------
  subroutine Iapbp_nn(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)

    Integer :: i,j, ip,jp, imin,imax, jmin,jmax
    Real(8) :: c

    call TimerStart('convol: I(a.;b.) n n')
    a(1:ns,1:ks)=0.d0
    do jp=1,ks+ks-1
       jmin=max(1,1+ks-jp)
       jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax
          c=0.d0
          do ip=1,ks+ks-1
             imin=max(1,1+ks-ip)
             imax=min(ns,ns+ks-ip)
             do i =imin,imax
                c=c+d(i,ip)*rkb(i,j,ip,jp)
             end do
          end do
          a(j,jp)=c
       end do
    end do
    call TimerStop('convol: I(a.;b.) n n')
  end subroutine Iapbp_nn
  ! --------------------------------------------------------------------------------
  subroutine Iapbp_ll(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)

    Integer :: i,j, ip,jp
    Real(8) :: c

    call TimerStart('convol: I(a.;b.) l l')
    a(1:ns,1:ks)=0.d0
    do jp=1,ks
       do j=ks+1-jp,ns
          c=0.d0
          do ip=1,ks
             do i=ks+1-ip,ns
                c=c+d(i,ip)*rkb(i,j,ip,jp)
             end do
          end do
          a(j,jp)=c
       end do
    end do
    call TimerStop('convol: I(a.;b.) l l')
  end subroutine Iapbp_ll
  ! --------------------------------------------------------------------------------
  subroutine Iapbp_ln(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)

    Integer :: i,j, ip,jp, jmin,jmax
    Real(8) :: c

    call TimerStart('convol: I(a.;b.) l n')
    a(1:ns,1:ks)=0.d0
    do jp=1,ks+ks-1
       jmin=max(1,1+ks-jp)
       jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax
          c=0.d0
          do ip=1,ks
             do i=ks+1-ip,ns
                c=c+d(i,ip)*rkb(i,j,ip,jp)
             end do
          end do
          a(j,jp)=c
       end do
    end do
    call TimerStop('convol: I(a.;b.) l n')
  end subroutine Iapbp_ln
  ! --------------------------------------------------------------------------------
  subroutine Iapbp_nl(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)

    Integer :: i,j, ip,jp, imin,imax
    Real(8) :: c

    call TimerStart('convol: I(a.;b.) n l')
    a(1:ns,1:ks)=0.d0
    do jp=1,ks
       do j=ks+1-jp,ns
          c=0.d0
          do ip=1,ks+ks-1
             imin=max(1,1+ks-ip)
             imax=min(ns,ns+ks-ip)
             do i =imin,imax
                c=c+d(i,ip)*rkb(i,j,ip,jp)
             end do
          end do
          a(j,jp)=c
       end do
    end do
    call TimerStop('convol: I(a.;b.) n l')
  end subroutine Iapbp_nl

  ! ----------------------------------------------------------------------------------------------------
  ! I(. a ; b .)
  ! ----------------------------------------------------------------------------------------------------
  subroutine Ipabp_ss(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    ! Use BLAS
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)

    Integer :: ip,jp
    ! Integer :: num_blas_threads
    Real(8) :: c

    ! call TimerStart('convol: I(.a;b.) s s, new')

    ! a(1:ns,1:ns) = 0.d0

    ! ! ! num_blas_threads = blas_get_num_threads()
    ! ! ! call blas_set_num_threads(1)
    ! ! ! !$omp parallel do collapse(2) default(none) private(ip,jp) shared(ks,ns,rkb,d24_rk)
    ! ! do ip=1,ks
    ! !    do jp=1,ks
    ! !       a( 1:ns-ks+1,1+ks-1:ns) = a( 1:ns-ks+1,1+ks-1:ns) + matmul(rkb(1:ns-ks+1,1:ns-ks+1,ip,jp), d(1+ks-1:ns, 1:ns-ks+1))

    ! !       if(ip.gt.1) then
    ! !          a(1+ks-1:ns,1+ks-1:ns) = a(1+ks-1:ns,1+ks-1:ns) + matmul(rkb(1:ns-ks+1,1:ns-ks+1,ip,jp), d( 1:ns-ks+1, 1:ns-ks+1))
    ! !       end if

    ! !       if(jp.gt.1) then
    ! !          a( 1:ns-ks+1, 1:ns-ks+1) = a( 1:ns-ks+1, 1:ns-ks+1) + matmul(rkb(1:ns-ks+1,1:ns-ks+1,ip,jp), d(1+ks-1:ns,1+ks-1:ns))
    ! !       end if

    ! !       if(ip.gt.1.and.jp.gt.1) then
    ! !          a(1+ks-1:ns, 1:ns-ks+1) = a(1+ks-1:ns, 1:ns-ks+1) + matmul(rkb(1:ns-ks+1,1:ns-ks+1,ip,jp), d( 1:ns-ks+1,1+ks-1:ns))
    ! !       end if
    ! !    end do
    ! ! end do
    ! ! ! !$omp end parallel do
    ! ! ! call blas_set_num_threads(num_blas_threads)
    ! call TimerStop('convol: I(.a;b.) s s, new')

#ifdef DEBUG_SPEEDUPS
    call Ipabp_ss_old(ns, ks, a, d)
#endif
  end subroutine Ipabp_ss


#ifdef DEBUG_SPEEDUPS
  subroutine Ipabp_ss_old(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in) :: d(ns,*)
    Real(8), intent(inout) :: a(ns,*)

    Integer :: i,j, ip,jp, ii,jj
    Real(8) :: c

    Real(8), allocatable :: a_ref(:,:)
    Real(8) :: discrepancy
    Real(8), parameter :: tolerance = sqrt(epsilon(1.d0))

    call TimerStart('convol: I(.a;b.) s s, old')
    allocate(a_ref(1:ns,1:ns))
    a_ref(1:ns,1:ns) = 0.d0
    do jp=1,ks
       do j =1,ns-jp+1
          jj=j+jp-1
          do ip=1,ks
             do i =1,ns-ip+1
                ii=i+ip-1
                c = rkb(i,j,ip,jp)
                a_ref( i,jj) = a_ref( i,jj) + c*d(ii, j)
                if(ip.gt.1)              a_ref(ii,jj) = a_ref(ii,jj) + c*d( i, j)
                if(jp.gt.1)              a_ref( i, j) = a_ref( i, j) + c*d(ii,jj)
                if(ip.gt.1.and.jp.gt.1)  a_ref(ii, j) = a_ref(ii, j) + c*d( i,jj)
             end do
          end do
       end do
    end do
    call TimerStop('convol: I(.a;b.) s s, old')
    ! discrepancy = sum(abs(a(1:ns,1:ns) - a_ref(1:ns,1:ns)))
    ! write(*,'("Discrepancy: ",e26.16,", tolerance: ",e26.16)') discrepancy, tolerance
    ! ! if(discrepancy > tolerance) error stop "Discrepancy exceeds tolerance"
    ! if(discrepancy > tolerance) then
    !    a(1:ns,1:ns) = a_ref(1:ns,1:ns)
    ! end if
    a(1:ns,1:ns) = a_ref(1:ns,1:ns)
    deallocate(a_ref)
  end subroutine Ipabp_ss_old
#endif
  ! --------------------------------------------------------------------------------
  subroutine Ipabp_ll(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in) :: d(ns,*)
    Real(8), intent(inout) :: a(ns,*)

    Integer :: i,j, ip,jp, ii,jj
    Real(8) :: c

    call TimerStart('convol: I(.a;b.) l l')
    do jp=1,ks
       do j=ks+1-jp,ns
          jj=j+jp-ks
          do ip=1,ks
             do i=ks+1-ip,ns
                ii=i+ip-ks
                c = rkb(i,j,ip,jp)
                a( i,jj) = a( i,jj) + c*d(ii, j)
                if(ip.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
                if(jp.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
                if(ip.ne.ks.and.jp.ne.ks) a(ii, j) = a(ii, j) + c*d( i,jj)
             end do
          end do
       end do
    end do
    call TimerStop('convol: I(.a;b.) l l')
  end subroutine Ipabp_ll
  ! --------------------------------------------------------------------------------
  subroutine Ipabp_nn(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in) :: d(ns,*)
    Real(8), intent(inout) :: a(ns,*)

    Integer :: i,j, ip,jp, ii,jj, imin,imax, jmin,jmax
    Real(8) :: c

    call TimerStart('convol: I(.a;b.) n n')
    do jp=1,ks+ks-1
       jmin=max(1,1+ks-jp)
       jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax
          jj=j+jp-ks
          do ip=1,ks+ks-1
             imin=max(1,1+ks-ip)
             imax=min(ns,ns+ks-ip)
             do i =imin,imax
                ii=i+ip-ks
                c = rkb(i,j,ip,jp)
                a( i,jj) = a( i,jj) + c*d(ii, j)
             end do
          end do
       end do
    end do
    call TimerStop('convol: I(.a;b.) n n')
  end subroutine Ipabp_nn
  ! --------------------------------------------------------------------------------
  subroutine Ipabp_nl(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in) :: d(ns,*)
    Real(8), intent(inout) :: a(ns,*)

    Integer :: i,j, ip,jp, ii,jj, imin,imax
    Real(8) :: c

    call TimerStart('convol: I(.a;b.) n l')
    do jp=1,ks
       do j=ks+1-jp,ns
          jj=j+jp-ks
          do ip=1,ks+ks-1
             imin=max(1,1+ks-ip)
             imax=min(ns,ns+ks-ip)
             do i =imin,imax
                ii=i+ip-ks
                c = rkb(i,j,ip,jp)
                a( i,jj) = a( i,jj) + c*d(ii, j)
                if(jp.ne.ks)             a( i, j) = a( i, j) + c*d(ii,jj)
             end do
          end do
       end do
    end do
    call TimerStop('convol: I(.a;b.) n l')
  end subroutine Ipabp_nl
  ! --------------------------------------------------------------------------------
  subroutine Ipabp_ln(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in) :: d(ns,*)
    Real(8), intent(inout) :: a(ns,*)

    Integer :: i,j, ip,jp, ii,jj, jmin,jmax
    Real(8) :: c

    call TimerStart('convol: I(.a;b.) l n')
    do jp=1,ks+ks-1
       jmin=max(1,1+ks-jp)
       jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax
          jj=j+jp-ks
          do ip=1,ks
             do i=ks+1-ip,ns
                ii=i+ip-ks
                c = rkb(i,j,ip,jp)
                a( i,jj) = a( i,jj) + c*d(ii, j)
                if(ip.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
             end do
          end do
       end do
    end do
    call TimerStop('convol: I(.a;b.) l n')
  end subroutine Ipabp_ln

  ! ----------------------------------------------------------------------------------------------------
  ! I(a . ; . b)
  ! ----------------------------------------------------------------------------------------------------
  subroutine Iappb_ss(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    ! Use BLAS
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in ) :: d(ns,*)
    Real(8), intent(out) :: a(ns,*)

    Integer :: ip,jp
    ! Integer :: num_blas_threads
    Real(8) :: c

    ! call TimerStart('convol: I(a.;.b) s s, new')

    ! a(1:ns,1:ns) = 0.d0
    ! call TimerStop('convol: I(a.;.b) s s, new')

#ifdef DEBUG_SPEEDUPS
    call Iappb_ss_old(ns, ks, a, d)
#endif
  end subroutine Iappb_ss


#ifdef DEBUG_SPEEDUPS
  subroutine Iappb_ss_old(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in) :: d(ns,*)
    Real(8), intent(inout) :: a(ns,*)

    Integer :: i,j, ip,jp, ii,jj
    Real(8) :: c

    Real(8), allocatable :: a_ref(:,:)
    Real(8) :: discrepancy
    Real(8), parameter :: tolerance = sqrt(epsilon(1.d0))

    call TimerStart('convol: I(a.;.b) s s, old')
    allocate(a_ref(1:ns,1:ns))
    a_ref(1:ns,1:ns) = 0.d0
    do jp=1,ks
       do j =1,ns-jp+1
          jj=j+jp-1
          do ip=1,ks
             do i =1,ns-ip+1
                ii=i+ip-1
                c = rkb(i,j,ip,jp)
                a_ref(ii, j) = a_ref(ii, j) + c*d( i,jj)
                if(ip.gt.1)              a_ref( i, j) = a_ref( i, j) + c*d(ii,jj)
                if(jp.gt.1)              a_ref(ii,jj) = a_ref(ii,jj) + c*d( i, j)
                if(ip.gt.1.and.jp.gt.1)  a_ref( i,jj) = a_ref( i,jj) + c*d(ii, j)
             end do
          end do
       end do
    end do
    call TimerStop('convol: I(a.;.b) s s, old')
    ! discrepancy = sum(abs(a(1:ns,1:ns) - a_ref(1:ns,1:ns)))
    ! write(*,'("Discrepancy: ",e26.16,", tolerance: ",e26.16)') discrepancy, tolerance
    ! ! if(discrepancy > tolerance) error stop "Discrepancy exceeds tolerance"
    ! if(discrepancy > tolerance) then
    !    a(1:ns,1:ns) = a_ref(1:ns,1:ns)
    ! end if
    a(1:ns,1:ns) = a_ref(1:ns,1:ns)
    deallocate(a_ref)
  end subroutine Iappb_ss_old
#endif
  ! --------------------------------------------------------------------------------
  subroutine Iappb_ll(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in) :: d(ns,*)
    Real(8), intent(inout) :: a(ns,*)

    Integer :: i,j, ip,jp, ii,jj
    Real(8) :: c

    call TimerStart('convol: I(a.;.b) l l')
    do jp=1,ks
       do j=ks+1-jp,ns
          jj=j+jp-ks
          do ip=1,ks
             do i=ks+1-ip,ns
                ii=i+ip-ks
                c = rkb(i,j,ip,jp)
                a(ii, j) = a(ii, j) + c*d( i,jj)
                if(ip.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
                if(jp.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
                if(ip.ne.ks.and.jp.ne.ks) a( i,jj) = a( i,jj) + c*d(ii, j)
             end do
          end do
       end do
    end do
    call TimerStop('convol: I(a.;.b) l l')
  end subroutine Iappb_ll
  ! --------------------------------------------------------------------------------
  subroutine Iappb_nn(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in) :: d(ns,*)
    Real(8), intent(inout) :: a(ns,*)

    Integer :: i,j, ip,jp, ii,jj, imin,imax, jmin,jmax
    Real(8) :: c

    call TimerStart('convol: I(a.;.b) n n')
    do jp=1,ks+ks-1
       jmin=max(1,1+ks-jp)
       jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax
          jj=j+jp-ks
          do ip=1,ks+ks-1
             imin=max(1,1+ks-ip)
             imax=min(ns,ns+ks-ip)
             do i =imin,imax
                ii=i+ip-ks
                c = rkb(i,j,ip,jp)
                a(ii, j) = a(ii, j) + c*d( i,jj)
             end do
          end do
       end do
    end do
    call TimerStop('convol: I(a.;.b) n n')
  end subroutine Iappb_nn
  ! --------------------------------------------------------------------------------
  subroutine Iappb_nl(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in) :: d(ns,*)
    Real(8), intent(inout) :: a(ns,*)

    Integer :: i,j, ip,jp, ii,jj, imin,imax
    Real(8) :: c

    call TimerStart('convol: I(a.;.b) n l')
    do jp=1,ks
       do j=ks+1-jp,ns
          jj=j+jp-ks
          do ip=1,ks+ks-1
             imin=max(1,1+ks-ip)
             imax=min(ns,ns+ks-ip)
             do i =imin,imax
                ii=i+ip-ks
                c = rkb(i,j,ip,jp)
                a(ii, j) = a(ii, j) + c*d( i,jj)
                if(jp.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
             end do
          end do
       end do
    end do
    call TimerStop('convol: I(a.;.b) n l')
  end subroutine Iappb_nl
  ! --------------------------------------------------------------------------------
  subroutine Iappb_ln(ns, ks, a, d)
    Use DBS_integrals, only: rkb
    Use Timer

    Implicit none

    Integer, intent(in) :: ns,ks
    Real(8), intent(in) :: d(ns,*)
    Real(8), intent(inout) :: a(ns,*)

    Integer :: i,j, ip,jp, ii,jj, jmin,jmax
    Real(8) :: c

    call TimerStart('convol: I(a.;.b) l n')
    do jp=1,ks+ks-1
       jmin=max(1,1+ks-jp)
       jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax
          jj=j+jp-ks
          do ip=1,ks
             do i=ks+1-ip,ns
                ii=i+ip-ks
                c = rkb(i,j,ip,jp)
                a(ii, j) = a(ii, j) + c*d( i,jj)
                if(ip.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
             end do
          end do
       end do
    end do
    call TimerStop('convol: I(a.;.b) l n')
  end subroutine Iappb_ln
End module DBS_convolutions
