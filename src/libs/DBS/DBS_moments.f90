!====================================================================
    Module DBS_moments
!====================================================================
!   contains moments defining as   <B_i|r^k|B_j>   over an interval.
!   These moments are used for calculation of two-electron integrals
!   according to the cell algorithm.
!--------------------------------------------------------------------
    Implicit none
    Integer :: kmk = -100, kk = 0
    Character(3) :: mtype = 'aaa'
    Real(8), allocatable :: rkd(:,:,:)
    Real(8), allocatable :: rkd1(:,:), rkd2(:,:), rkd3(:,:), rkd4(:,:)
    Real(8) :: memory_DBS_moments

    End Module DBS_moments

!====================================================================
    Subroutine alloc_DBS_moments
!====================================================================
!   allocate arrays in module DBS_moments and initialize the flags
!   mtype and kmk
!-------------------------------------------------------------------
    Use DBS_grid, only: nv,ks
    Use DBS_moments

    if(allocated(rkd)) Deallocate(rkd,rkd1,rkd2,rkd3,rkd4)
    kk = ks*ks
    Allocate(rkd(kk,kk,nv), rkd1(kk,nv),rkd2(kk,nv), &
                            rkd3(kk,nv),rkd4(kk,nv))
    rkd = 0.d0; rkd1 = 0.d0; rkd2 = 0.d0; rkd3 = 0.d0; rkd4 = 0.d0

    mtype='bbb';  kmk=-100

    memory_DBS_moments = (kk*kk*nv + 4*kk*nv)*8.d0 / (1024*1024)

    End Subroutine alloc_DBS_moments


!====================================================================
    Subroutine dealloc_DBS_moments
!====================================================================
!   deallocate array in module DBS_moments and initialize the flags
!   mtype and kmk
!--------------------------------------------------------------------
    Use DBS_moments
    if(allocated(rkd)) Deallocate(rkd, rkd1,rkd2,rkd3,rkd4)
    mtype='aaa';  kmk=-100
    End Subroutine dealloc_DBS_moments


!=====================================================================
    Subroutine  moments (k,rkm,sym,dir)
!=====================================================================
!   Computes moments defining as <B_i|r^k|B_j> over an interval
!---------------------------------------------------------------------
!   on entry
!   --------
!      k      the power of moments
!      sym    's' or 'n' - symmetrical or non-symetrical case
!      div    'b' or 'd' - B_j or B'_j
!   on exit
!   -------
!      rkm    array of moments over all intervals.
!---------------------------------------------------------------------
      Use DBS_grid

      Implicit none
      Integer, intent(in) :: k
      Real(8), intent(out) :: rkm(ks*ks,nv)
      Character, intent(in) :: sym, dir
      ! local variables
      Integer :: iv, jk
      Real(8) :: hp
      Real(8) :: rv(ks*ks)

      jk = ks*(ks+1)/2; if(sym.eq.'n') jk = ks*ks

      if(me.eq.0) then     ! .. there is no exponential grid

       Do iv = 1,nv
        Call moment(iv,k,rv,sym,dir);  rkm(1:jk,iv) = rv(1:jk)
       End do

      else                 ! .. case of exponential grid

       ! .. the first non-exponential region
       DO iv = 1,ml+ks-1
        Call moment(iv,k,rv,sym,dir);  rkm(1:jk,iv) = rv(1:jk)
       End DO

       ! .. the exponential region -> using scaling law
       hp = (1.d0+he)**k; if(dir.eq.'b') hp = hp*(1.d0+he)
       DO iv=ml+ks,ml+me-ks+2
        rkm(1:jk,iv) = rkm(1:jk,iv-1)*hp
       End DO

       ! .. the last non-exponential region
       DO iv=ml+me-ks+3,nv
        Call moment(iv,k,rv,sym,dir);  rkm(1:jk,iv) = rv(1:jk)
       End DO

      End if

      End Subroutine moments


!=====================================================================
    Subroutine  moment (iv,k,rv,sym,dir)
!=====================================================================
!   Computes moment defining as <B_i|r^k|B_j> for given interval
!---------------------------------------------------------------------
!   on entry
!   --------
!      k      the power of moments
!      sym    's' or 'n' - symmetrical or non-symetrical case
!      div    'b' or 'd' - B_j or B'_j
!      iv     index of interval
!   on exit
!   -------
!      rv     array of moments for interval 'iv'
!---------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: k,iv
      Real(8), intent(out) :: rv(ks*ks)
      Character, intent(in) :: sym, dir
      ! local variables
      Integer :: i, j, ii, jj
      Real(8) :: bi(ks,ks),bj(ks,ks)

      gw = grw(iv,:)
      if( k > 0 ) gw = grw(iv,:) * gr (iv,:)**(+k)
      if( k < 0 ) gw = grw(iv,:) * grm(iv,:)**(-k)

      bi(:,:) = bsp(iv,:,:)
      if(dir.eq.'b') then
       Do j = 1,ks;  bj(:,j) = bsp(iv,:,j)*gw(:); End do
      else
       Do j = 1,ks;  bj(:,j) = bsq(iv,:,j)*gw(:); End do
      End if

      ii = 0
      Do i=1,ks; jj = 1; if(sym.eq.'s') jj = i
       Do j=jj,ks
        ii = ii + 1;  rv(ii) = SUM(bi(:,i)*bj(:,j))
       End do
      End do

      End Subroutine moment


!======================================================================
      Subroutine  moments_pp(k,kk,nv,rkm)
!======================================================================
!     Computes moments between P-functions
!     (cheking the scaling properties in future if lucky)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k,kk,nv
      Real(8), intent(out) :: rkm(kk,nv)
      Integer :: iv

      rkm = 0.d0
      Do iv = 1,nv; Call moment_pp(iv,k,rkm(1,iv)); End do

      End Subroutine moments_pp


!======================================================================
      Subroutine moment_pp(iv,k,rv)
!======================================================================
!     Computes moment between P-functions for interval "iv"
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: iv,k
      Real(8), intent(out) :: rv(*)
      Integer :: i, j, ii
      Real(8) :: b1(ks,ks),b2(ks,ks)

      gw = grw(iv,1:ks)
      if( k > 0 ) gw = grw(iv,1:ks) * gr (iv,1:ks)**(+k)
      if( k < 0 ) gw = grw(iv,1:ks) * grm(iv,1:ks)**(-k)

      b1(1:ks,1:ksp) = pbsp(iv,1:ks,1:ksp)
      Do j=1,ksp; b2(1:ks,j)=b1(1:ks,j)*gw(1:ks); End do

      ii = 0
      Do i=1,ksp; Do j=i,ksp
        ii=ii+1; rv(ii)=SUM(b1(1:ks,i)*b2(1:ks,j))
      End do; End do

      End Subroutine moment_pp


!======================================================================
      Subroutine  moments_qq(k,kk,nv,rkm)
!======================================================================
!     Computes moments between Q-functions
!     (cheking the scaling properties in future if lucky)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k,kk,nv
      Real(8), intent(out) :: rkm(kk,nv)
      Integer :: iv

      rkm = 0.d0
      Do iv = 1,nv; Call moment_qq(iv,k,rkm(1,iv)); End do

      End Subroutine moments_qq


!======================================================================
      Subroutine moment_qq(iv,k,rv)
!======================================================================
!     Computes moment between Q-functions for interval "iv"
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: iv,k
      Real(8), intent(out) :: rv(*)
      Integer :: i, j, ii
      Real(8) :: b1(ks,ks),b2(ks,ks)

      gw = grw(iv,1:ks)
      if( k > 0 ) gw = grw(iv,1:ks) * gr (iv,1:ks)**(+k)
      if( k < 0 ) gw = grw(iv,1:ks) * grm(iv,1:ks)**(-k)

      b1(1:ks,1:ksq) = qbsp(iv,1:ks,1:ksq)
      Do j=1,ksq; b2(1:ks,j)=b1(1:ks,j)*gw(1:ks); End do

      ii = 0
      Do i=1,ksq; Do j=i,ksq
        ii=ii+1; rv(ii)=SUM(b1(1:ks,i)*b2(1:ks,j))
      End do; End do

      End Subroutine moment_qq


!======================================================================
      Subroutine  moments_pq(k,kk,nv,rkm)
!======================================================================
!     Computes moments between P and Q-functions
!     (cheking the scaling properties in future if lucky)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k,kk,nv
      Real(8), intent(out) :: rkm(kk,nv)
      Integer :: iv

      rkm = 0.d0
      Do iv = 1,nv; Call moment_pq(iv,k,rkm(1,iv)); End do

      End Subroutine moments_pq


!======================================================================
      Subroutine moment_pq(iv,k,rv)
!======================================================================
!     Computes moment between P- and Q-functions for interval "iv"
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: iv,k
      Real(8), intent(out) :: rv(*)
      Integer :: i, j, ii
      Real(8) :: b1(ks,ks),b2(ks,ks)

      gw = grw(iv,1:ks)
      if( k > 0 ) gw = grw(iv,1:ks) * gr (iv,1:ks)**(+k)
      if( k < 0 ) gw = grw(iv,1:ks) * grm(iv,1:ks)**(-k)

      b1(1:ks,1:ksp) = pbsp(iv,1:ks,1:ksp)
      Do j=1,ksq; b2(1:ks,j)=qbsp(iv,1:ks,j)*gw(1:ks); End do

      ii = 0
      Do i=1,ksp; Do j=1,ksq
        ii=ii+1; rv(ii)=SUM(b1(1:ks,i)*b2(1:ks,j))
      End do; End do

      End Subroutine moment_pq


!======================================================================
      Subroutine  moments_qp(k,kk,nv,rkm)
!======================================================================
!     Computes moments between Q- and P-functions
!     (cheking the scaling properties in future if lucky)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k,kk,nv
      Real(8), intent(out) :: rkm(kk,nv)
      Integer :: iv

      rkm = 0.d0
      Do iv = 1,nv; Call moment_qp(iv,k,rkm(1,iv)); End do

      End Subroutine moments_qp


!======================================================================
      Subroutine moment_qp(iv,k,rv)
!======================================================================
!     Computes moment between Q- and P-functions for interval "iv"
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: iv,k
      Real(8), intent(out) :: rv(*)
      Integer :: i, j, ii
      Real(8) :: b1(ks,ks),b2(ks,ks)

      gw = grw(iv,1:ks)
      if( k > 0 ) gw = grw(iv,1:ks) * gr (iv,1:ks)**(+k)
      if( k < 0 ) gw = grw(iv,1:ks) * grm(iv,1:ks)**(-k)

      b1(1:ks,1:ksq) = qbsp(iv,1:ks,1:ksq)
      Do j=1,ksp; b2(1:ks,j)=pbsp(iv,1:ks,j)*gw(1:ks); End do

      ii = 0
      Do i=1,ksq; Do j=1,ksp
        ii=ii+1; rv(ii)=SUM(b1(1:ks,i)*b2(1:ks,j))
      End do; End do

      End Subroutine moment_qp

