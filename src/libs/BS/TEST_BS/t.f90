!======================================================================
    PROGRAM test_bs
!======================================================================
!
!   This program sets up some basic arrays in the spline basis,
!   determines a set of hydrogenic orbitals,
!   and evaluates a series of radial integrals
!
!   SUBROUTINE called:
!       define_grid
!       define_spline
!
!   CONTAINS:
!       define_orbitals
!       test_integrals
!
!----------------------------------------------------------------------

    USE spline_atomic

    IMPLICIT NONE

    ! .. sets up positions of grid points: t

    CALL define_grid (z)

    ! .. initializes the values of the spline and its derivatives
    ! .. and evaluates the spline arrays (operators in spline basis)
    ! .. which are defined in the MODULES spline_grid and spline_galerkin

    CALL define_spline


    ! .. defines the set of hydrogenic orbitals

    Call define_orbitals

    fine = 1.d0                  ! test assignment

    CALL test_integrals

    Call test_quadr(z)

    Call test_hl(z)

    Call test_azl(z)

    Call test_vc(z)


    END PROGRAM test_bs


!========================================================================
      SUBROUTINE define_orbitals
!========================================================================

! ... define the set of hydrogenic orbitals for n<=5, l<=4 (15 at all)

      USE spline_param
      USE spline_atomic
      USE spline_orbitals

      IMPLICIT NONE

      ! .. local variables

      INTEGER :: j,nj,lj
      CHARACTER (LEN=4) :: ELF4
      REAL(KIND=8), DIMENSION(ns) :: coef

      Call allocate_bsorb(15)
      nbf = mbf

      j=0
      Do nj=1,5
        Do lj=0,nj-1
         z = 2.0; if(nj.gt.1) z=1
          CALL bhwf(nj,lj,z,coef)
          j = j+1
          pbs(:,j) = coef(:)
          nbs(j) = nj
          lbs(j) = lj
          mbs(j) = ns
          kbs(j) = 0
          Ebs(j) = ELF4(nj,lj,0)
        End do
      End do

      END SUBROUTINE define_orbitals



!========================================================================
      SUBROUTINE test_integrals
!========================================================================
!
!     tests the given set of radial integrals on the basis of hydrogenic
!     orbitals. The integrals are given in file 'int_inp', the results are
!     output in file 'int_out'.
!
!------------------------------------------------------------------------

      USE spline_param
      USE spline_orbitals
      USE spline_atomic
      USE spline_integrals

      IMPLICIT REAL(8) (A-H,O-Z)

      Character INT*1, AAA*7, AS*80
      Character( 4) :: EL1, EL2, EL3, EL4

      INTEGER(4) :: i,k, in,iout

      INTEGER(4) :: n1,n2,n3,n4
      INTEGER(4) :: l1,l2,l3,l4
      INTEGER(4) :: k1,k2,k3,k4
      INTEGER(4) :: i1,i2,i3,i4
      INTEGER(4) :: ie1,ie2,ie3,ie4

      REAL(8) :: S, SS, SA,SB, S1, S2

      REAL(8) :: RK,RKc,RKy,RKd,RKe,RKv
      REAL(8) :: NK,NKc,NKy,NKd,NKe,NKv
      REAL(8) :: MK,MKc,MKy,MKd,MKe,MKv
      REAL(8) :: TK,TKc,TKy,TKd,TKe,TKv
      REAL(8) :: VK,VKc,VKy,VKd,VKe,VKv
      REAL(8) :: WK,WKc,WKy,WKd,WKe,WKv
      REAL(8) :: QK,QKc,QKy,QKd,QKe,QKv

      REAL(8) :: T1,T2

      INTEGER(4), EXTERNAL :: Ifind_bsorb
      REAL(8), EXTERNAL :: RRTC

! ... input ile

      in=10
      Open(in,file='int_inp')

! ... output file

      iout=11
      Open(iout,file='int_out')

! ... loop over input strings

    1 read(in,'(a)',end=2) AS

      if(AS(1:1).eq.'*') go to 2          ! end of data

      if(AS(1:1).eq.'#') then             ! timing block
!----------------------------------------------------------------------

      INT = AS(2:2)

! ... timing for mrk_diff

      T1 = RRTC()
      krk=-100
      Do k = 0,4
       if(INT.eq.'R') Call mrk_diff(k)
       if(INT.eq.'N') Call mnk_diff(k)
       if(INT.eq.'M') Call mmk_diff(k)
       if(INT.eq.'T') Call mtk_diff(k)
       if(INT.eq.'V') Call mvk_diff(k)
      End do
      T_diff = (RRTC() - T1)/5

! ... timing for mrk_cell

      T1 = RRTC()
      krk=-100
      Do k = 0,4
       if(INT.eq.'R') Call mrk_cell(k)
       if(INT.eq.'N') Call mnk_cell(k)
       if(INT.eq.'M') Call mmk_cell(k)
       if(INT.eq.'T') Call mtk_cell(k)
       if(INT.eq.'V') Call mvk_cell(k)
       if(INT.eq.'W') Call mwk_cell(k)
       if(INT.eq.'Q') Call mqk_cell(k)
      End do
      T_cell = (RRTC() - T1)/5

! ... timing for rk_moment

      T1 = RRTC()
      Call rk_moments(5)
      Do k = 0,4
       if(INT.eq.'R') Call rk_moments(k)
       if(INT.eq.'N') Call nk_moments(k)
       if(INT.eq.'M') Call mk_moments(k)
       if(INT.eq.'T') Call tk_moments(k)
       if(INT.eq.'V') Call vk_moments(k)
       if(INT.eq.'V') Call wk_moments(k)
       if(INT.eq.'Q') Call qk_moments(k)
      End do
      T_moments = (RRTC() - T1)/5

! ... timing for RKy

      ie1 = Ifind_bsorb(2,1,0); ie2 = Ifind_bsorb(3,2,0)
      T1 = RRTC()
      Do i=1,1000
       if(INT.eq.'R')  S = RKy (ie1,ie1,ie1,ie1,0)
       if(INT.eq.'R')  S = RKy (ie2,ie2,ie2,ie2,2)
       if(INT.eq.'N')  S = NKy (ie1,ie1,ie1,ie1,0)
       if(INT.eq.'N')  S = NKy (ie2,ie2,ie2,ie2,2)
       if(INT.eq.'M')  S = MKy (ie1,ie1,ie1,ie1,0)
       if(INT.eq.'M')  S = MKy (ie2,ie2,ie2,ie2,2)
       if(INT.eq.'T')  S = TKy (ie1,ie1,ie1,ie1,0)
       if(INT.eq.'T')  S = TKy (ie2,ie2,ie2,ie2,2)
       if(INT.eq.'V')  S = VKy (ie1,ie1,ie1,ie1,0)
       if(INT.eq.'V')  S = VKy (ie2,ie2,ie2,ie2,2)
       if(INT.eq.'W')  S = WKy (ie1,ie1,ie1,ie1,0)
       if(INT.eq.'W')  S = WKy (ie2,ie2,ie2,ie2,2)
       if(INT.eq.'Q')  S = QKy (ie1,ie1,ie1,ie1,0)
       if(INT.eq.'Q')  S = QKy (ie2,ie2,ie2,ie2,2)
      End do
      T_y = (RRTC()-T1)/2000

! ... timing for assembling B-spline integrals

      k = 0; Call mrk_cell(k); ie1 = Ifind_bsorb(3,0,0)

      T1 = RRTC()
      Do i=1,500
        if(INT.eq.'R') S = RK (ie1,ie1,ie1,ie1,k,'c')
        if(INT.eq.'N') S = NK (ie1,ie1,ie1,ie1,k,'c')
        if(INT.eq.'M') S = MK (ie1,ie1,ie1,ie1,k,'c')
        if(INT.eq.'T') S = TK (ie1,ie1,ie1,ie1,k,'c')
        if(INT.eq.'V') S = VK (ie1,ie1,ie1,ie1,k,'c')
        if(INT.eq.'W') S = WK (ie1,ie1,ie1,ie1,k,'c')
        if(INT.eq.'Q') S = QK (ie1,ie1,ie1,ie1,k,'c')
      End do
      T_a = (RRTC()-T1)/500

! ... timing for assembling the moments

      T1 = RRTC()
      Do i=1,500
        if(INT.eq.'R') S = RKc(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'N') S = NKc(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'M') S = MKc(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'T') S = TKc(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'V') S = VKc(ie1,ie1,ie1,ie1,k)
      End do
      T_m = (RRTC()-T1)/500

! ... timing for convolution routins

      T1 = RRTC()
      Do i=1,500
        if(INT.eq.'R') S = RKd(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'N') S = NKd(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'M') S = MKd(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'T') S = TKd(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'V') S = VKd(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'W') S = WKd(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'Q') S = QKd(ie1,ie1,ie1,ie1,k)
      End do
      T_d = (RRTC()-T1)/500

      T1 = RRTC()
      Do i=1,500
        if(INT.eq.'R') S = RKe(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'N') S = NKe(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'M') S = MKe(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'T') S = TKe(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'V') S = VKe(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'W') S = WKe(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'Q') S = QKe(ie1,ie1,ie1,ie1,k)
      End do
      T_e = (RRTC()-T1)/500

      T1 = RRTC()
      Do i=1,500
        if(INT.eq.'R') S = RKv(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'N') S = NKv(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'M') S = MKv(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'T') S = TKv(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'V') S = VKv(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'W') S = WKv(ie1,ie1,ie1,ie1,k)
        if(INT.eq.'Q') S = QKv(ie1,ie1,ie1,ie1,k)
      End do
      T_v = (RRTC()-T1)/500


      write(iout,'(/a,a,a/)') &
       ' timing for ', INT, ' integrals (msec) :'
      write(iout,'(a,f10.3)') &
       ' generating the B-spline integrals  (_DIFF): ', T_diff
      write(iout,'(a,f10.3)') &
       ' generating the B-spline integrals  (_CELL): ', T_cell
      write(iout,'(a,f10.3)') &
       ' generating the B-spline moments (_MOMENTS): ', T_moments

      write(iout,'(a,2f10.3)') &
       ' differential equation method:   ',T_y, T_y/T_y
      write(iout,'(a,2f10.3)') &
       ' assembling B-spline integrals:  ',T_a, T_a/T_y
      write(iout,'(a,2f10.3)') &
       ' assembling of moments:          ',T_m, T_m/T_y
      write(iout,'(a,2f10.3)') &
       ' direct convolution:             ',T_d, T_d/T_y
      write(iout,'(a,2f10.3)') &
       ' exhange convolution:            ',T_e, T_e/T_y
      write(iout,'(a,2f10.3)') &
       ' vector convolution:             ',T_v, T_v/T_y

      go to 1
      end if                ! end of timing mode (#)


!----------------------------------------------------------------------
!
!         ... test of accuracy ...


      Read(AS,'(1x,i2,4(1x,a3))') k, el1, el2, el3, el4

      Call EL4_nlk(el1,n1,l1,k1);    ie1 = Ifind_bsorb(n1,l1,k1)
      Call EL4_nlk(el2,n2,l2,k2);    ie2 = Ifind_bsorb(n2,l2,k2)
      Call EL4_nlk(el3,n3,l3,k3);    ie3 = Ifind_bsorb(n3,l3,k3)
      Call EL4_nlk(el4,n4,l4,k2);    ie4 = Ifind_bsorb(n4,l4,k2)

      read(AS(22:80),*) S1,S2;  SB = S1/S2
      INT = AS(1:1);  SA = SB * Z**3; if(INT.eq.'R') SA = SB * Z

      write(iout,'(a21,d24.16)') AS(1:21), SA

      if(INT.eq.'R') S = RKY (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'N') S = NKY (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'M') S = MKY (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'T') S = TKY (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'V') S = VKY (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'W') S = WKY (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'Q') S = QKY (ie1,ie2,ie3,ie4,k)
      SS=S; if(SA.ne.0.d0) SS = (S-SA)/SA; AAA = INT// 'Ky'
      write(iout,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      write(*   ,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA

      if(INT.eq.'R') S = RKc (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'N') S = NKc (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'M') S = MKc (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'T') S = TKc (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'V') S = VKc (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'W') S = 0.d0
      if(INT.eq.'Q') S = 0.d0
      SS=S; if(SA.ne.0.d0) SS = (S-SA)/SA; AAA = INT// 'Kc'
      write(iout,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      write(*   ,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA

      krk = -100
      if(INT.eq.'R') S = RK  (ie1,ie2,ie3,ie4,k,'d')
      if(INT.eq.'N') S = NK  (ie1,ie2,ie3,ie4,k,'d')
      if(INT.eq.'M') S = MK  (ie1,ie2,ie3,ie4,k,'d')
      if(INT.eq.'T') S = TK  (ie1,ie2,ie3,ie4,k,'d')
      if(INT.eq.'V') S = VK  (ie1,ie2,ie3,ie4,k,'d')
      if(INT.eq.'W') S = 0.d0
      if(INT.eq.'Q') S = 0.d0
      SS=S; if(SA.ne.0.d0) SS = (S-SA)/SA; AAA = INT// 'K_diff'
      write(iout,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      write(*   ,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA

      krk = -100
      if(INT.eq.'R') S = RK  (ie1,ie2,ie3,ie4,k,'c')
      if(INT.eq.'N') S = NK  (ie1,ie2,ie3,ie4,k,'c')
      if(INT.eq.'M') S = MK  (ie1,ie2,ie3,ie4,k,'c')
      if(INT.eq.'T') S = TK  (ie1,ie2,ie3,ie4,k,'c')
      if(INT.eq.'V') S = VK  (ie1,ie2,ie3,ie4,k,'c')
      if(INT.eq.'W') S = WK  (ie1,ie2,ie3,ie4,k,'c')
      if(INT.eq.'Q') S = QK  (ie1,ie2,ie3,ie4,k,'c')
      SS=S; if(SA.ne.0.d0) SS = (S-SA)/SA; AAA = INT// 'K_cell'
      write(iout,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      write(*   ,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA

      if(INT.eq.'R') S = RKd (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'N') S = NKd (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'M') S = MKd (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'T') S = TKd (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'V') S = VKd (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'W') S = WKd (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'Q') S = QKd (ie1,ie2,ie3,ie4,k)
      SS=S; if(SA.ne.0.d0) SS = (S-SA)/SA; AAA = INT// 'Kd'
      write(iout,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      write(*   ,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA

      if(INT.eq.'R') S = RKe (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'N') S = NKe (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'M') S = MKe (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'T') S = TKe (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'V') S = VKe (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'W') S = WKe (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'Q') S = QKe (ie1,ie2,ie3,ie4,k)
      SS=S; if(SA.ne.0.d0) SS = (S-SA)/SA; AAA = INT// 'Ke'
      write(iout,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      write(*   ,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA

      if(INT.eq.'R') S = RKv (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'N') S = NKv (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'M') S = MKv (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'T') S = TKv (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'V') S = VKv (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'W') S = WKv (ie1,ie2,ie3,ie4,k)
      if(INT.eq.'Q') S = QKv (ie1,ie2,ie3,ie4,k)
      SS=S; if(SA.ne.0.d0) SS = (S-SA)/SA; AAA = INT// 'Kv'
      write(iout,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      write(*   ,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA

      if(INT.eq.'R') S = RKdn(ie1,ie2,ie3,ie4,k)
      SS=S; if(SA.ne.0.d0) SS = (S-SA)/SA; AAA = INT// 'Nd'
      write(iout,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      write(*   ,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      if(INT.eq.'R') S = RKen(ie1,ie2,ie3,ie4,k)
      SS=S; if(SA.ne.0.d0) SS = (S-SA)/SA; AAA = INT// 'Ne'
      write(iout,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      write(*   ,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      if(INT.eq.'R') S = RKvn(ie1,ie2,ie3,ie4,k)
      SS=S; if(SA.ne.0.d0) SS = (S-SA)/SA; AAA = INT// 'Nv'
      write(iout,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      write(*   ,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      if(INT.eq.'R') S = RK1n(ie1,ie2,ie3,ie4,k)
      SS=S; if(SA.ne.0.d0) SS = (S-SA)/SA; AAA = INT// 'N1'
      write(iout,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA
      write(*   ,'(21x,d24.16,d12.3,5x,a)') s,ss,AAA

      go to 1    ! go to new integral

    2 Continue   ! end of data

      END SUBROUTINE test_integrals




!========================================================================
      SUBROUTINE test_quadr(z)
!========================================================================
!
!     tests the accuracy of one-electron radial integrals on the basis
!     of hydrogenic orbitals. Results are output in file 'quadr_out'.
!
!------------------------------------------------------------------------

      USE spline_orbitals

      IMPLICIT NONE

      REAL(KIND=8), INTENT(in) :: Z
      REAL(KIND=8), EXTERNAL :: QUADR, GRAD
      INTEGER, EXTERNAL :: Ifind_bsorb
      INTEGER :: iout, i, j, icase
      REAL(KIND=8) :: a,ar1,ar2,ar3,am1,am2,am3,br1,br2,br3,bm1,bm2,bm3

      iout=2
      Open(iout,file='quadr_out')

!----------------------------------------------------------------------
!                                            accuracy of orthogonality:

      write(iout,'(/a/)')   '  accuracy of orthogonality'
      Do i=1,nbf
       Do j=i,nbf

       if(lbs(i).ne.lbs(j)) Cycle

       a = QUADR(i,j,0)
       if(i.eq.j) a = a - 1.d0

       write(iout,'(3x,a,3x,a,D12.3)') EBS(i),EBS(j),a

       End do
      End do

!-----------------------------------------------------------------------
!                                           accuracy of avarages <r^m>:

      write(iout,'(/a/)')   '  accuracy of <r^m>'
      write(iout,'(/4x,6(7x,a5)/)') &
          ' <r> ','<r^2>','<r^3>','<r-1>','<r-2>','<r-3>'

      Do i=1,nbf
       write(iout,*)
       Call HVA(nbs(i),lbs(i),z,ar1,ar2,ar3,am1,am2,am3)
       br1 = QUADR(i,i,1)-ar1
       br2 = QUADR(i,i,2)-ar2
       br3 = QUADR(i,i,3)-ar3
       bm1 = QUADR(i,i,-1)-am1
       bm2 = QUADR(i,i,-2)-am2
       bm3 = 0.0
       if(lbs(i).gt.0) bm3 = QUADR(i,i,-3)-am3
       if(ar1.ne.0.d0) br1=br1/ar1
       if(ar2.ne.0.d0) br2=br2/ar2
       if(ar3.ne.0.d0) br3=br3/ar3
       if(am1.ne.0.d0) bm1=bm1/am1
       if(am2.ne.0.d0) bm2=bm2/am2
       if(am3.ne.0.d0) bm3=bm3/am3
       write(iout,'(a4,1x,6D12.3)') ebs(i),ar1,ar2,ar3,am1,am2,am3
       write(iout,'(   5x,6D12.3)')        br1,br2,br3,bm1,bm2,bm3
      end do

!-----------------------------------------------------------------------
!                                     accuracy of transition amplitudes:

       write(iout,'(/a/)')  &
              '  accuracy of <i|r|j>:             QUADR  and  GRAD'

       Do icase = 1,3
        if(icase.eq.1)   i = Ifind_bsorb(1,0,0)         ! 1s - np
        if(icase.eq.2)   i = Ifind_bsorb(2,0,0)         ! 2s - np
        if(icase.eq.3)   i = Ifind_bsorb(2,1,0)         ! 2p - ns,nd

       Do j=1,nbf
        if(ABS(lbs(i)-lbs(j)).ne.1) Cycle
        if(nbs(j).le.nbs(i)) Cycle
        a = nbs(j)

        if(icase.eq.1) &
           ar1 = 16*a**3.5d0*(a-1)**(a-2.5d0)/(a+1)**(a+2.5d0)
        if(icase.eq.2) &
           ar1 = 2**8*a**3.5d0*sqrt(2*(a**2-1))*(a-2)**(a-3)/(a+2)**(a+3)
        if(icase.eq.3.and.lbs(j).eq.0) &
           ar1 = 2**8*a**4.5d0/sqrt(6.d0)*(a-2)**(a-3)/(a+2)**(a+3)
        if(icase.eq.3.and.lbs(j).eq.2) &
           ar1 = 2**10*a**4.5d0*sqrt((a**2-1)/6)* &
                 (a-2)**(a-3.5d0)/(a+2)**(a+3.5d0)

        ar1 =  ar1/z
        br1 =  (QUADR(i,j,1) - ar1)/ar1
        a   =  (1.d0/nbs(i)**2 - 1.d0/nbs(j)**2)/2.d0 * z**2
        bm1 =  GRAD(i,j) / a
        bm1 =  (bm1 - ar1)/ar1
        write(iout,'(a,a,a,a,a,3D12.3)') &
                    '<',EBS(i),'| r |',EBS(j),'> = ',ar1,br1,bm1
       End do
       End do

     END SUBROUTINE test_quadr





!========================================================================
      SUBROUTINE test_hl(z)
!========================================================================

! ... tests the accuracy of one-electron L integrals on the basis
! ... of hydrogenic orbitals. Results are output in file 'hl_out'.

      USE spline_orbitals

      IMPLICIT NONE

      REAL(KIND=8), INTENT(in) :: Z
      REAL(KIND=8), EXTERNAL :: BHL
      INTEGER :: iout, i, j
      REAL(KIND=8) :: a,b,c

      iout=2
      Open(iout,file='hl_out')

      write(iout,'(/a/)') '   accuracy and symmetry of HL:'

      Do i=1,nbf
       Do j=i,nbf

       if(lbs(i).ne.lbs(j)) Cycle

       a=0.d0
       if(i.eq.j) a = z*z/nbs(i)**2

       b=BHL(i,j)
       c=BHL(j,i)

!       if(a.ne.0.d0) b=(b-a)/a
!       if(a.ne.0.d0) c=(c-a)/a
!       write(iout,'(a,a,a,a,a,3D12.3)')  &
!                    '<',EBS(i),'| L |',EBS(j),'> =',a,b,c

       write(iout,'(a,a,a,a,a,2D20.8)')  &
                    '<',EBS(i),'| L |',EBS(j),'> =',b,c

       End do
      End do

      Close(iout)

    END SUBROUTINE test_hl



!========================================================================
      SUBROUTINE test_azl(z)
!========================================================================

! ... tests the accuracy of AZL subprogram on the basis
! ... of hydrogenic orbitals. Results are output in file 'azl_out'.

      USE spline_param
      USE spline_galerkin
      USE spline_grid
      USE spline_orbitals
      USE spline_hl

      IMPLICIT NONE

      REAL(KIND=8), INTENT(in) :: Z
      REAL(KIND=8), EXTERNAL :: HNORM, BVMV, BVALU2
      INTEGER :: iout, i,j, l,ns1, ii, jj
      REAL(KIND=8) :: a,b,ai,bi,aj,bj, b2s, b2ss, b3ss,fine
      REAL(KIND=8), DIMENSION(nv,ks,ks) :: dbiatx
      REAL(KIND=8), DIMENSION(nv) :: x,y,y1,y2
      REAL(KIND=8), DIMENSION(ns) :: bcoef, eval
      REAL(KIND=8), DIMENSION(ns,ns) :: hm,c
      REAL(8), External :: AZL


      INTEGER :: Lwork, info
      REAL(KIND=8), DIMENSION(3*ns) :: work
      Lwork = 3*ns

      iout=2
      Open(iout,file='azl_out')

      write(iout,'(/a/)') '   accuracy of the AZ coefficients: '

      Do i=1,nbf

        A = HNORM(nbs(i),lbs(i),z) * (2.*z/nbs(i))**(lbs(i) + 1)

        B = azl(z,h,ks,lbs(i)+1) * pbs(lbs(i)+2,i)

        write(iout,'(a,a,3D20.8)') EBS(i),':   ', a, b, (a-b)/a

      End do

      write(iout,'(/a/)') '   accuracy of ove-body Darwin terms: '

      Do i=1,nbf
       if(lbs(i).ne.0) Cycle
       ai = HNORM(nbs(i),lbs(i),z) * (2.*z/nbs(i))**(lbs(i) + 1)
       bi = azl(z,h,ks,lbs(i)+1) * pbs(lbs(i)+2,i)

       Do j=1,nbf
        if(lbs(j).ne.0) Cycle
        aj = HNORM(nbs(j),lbs(j),z) * (2.*z/nbs(j))**(lbs(j) + 1)
        bj = azl(z,h,ks,lbs(j)+1) * pbs(lbs(j)+2,j)

        a = ai * aj * z
        b = bi * bj * z

        write(iout,'(a4,a4,a,3D20.8)') &
                     EBS(i),EBS(j),': ', a, b, (a-b)/a

        hl = 0.d0
        b = azl(z,h,ks,lbs(i)+1)
        hl(2,ks) =  z*b*b
        b = BVMV (ns,ks,hl,'s',pbs(1,i),pbs(1,j))
        write(iout,'(a4,a4,a,3D20.8)') &
                     EBS(i),EBS(j),': ', a, b, (a-b)/a

       End do
      End do


      write(iout,'(/a/)') 'accuracy of ove-body Darwin terms: 2z^4/n^3'

      Do i=1,nbf

       if(lbs(i).ne.0) Cycle

       a = 2*z**4/nbs(i)**3

       bi = azl(z,h,ks,lbs(i)+1) * pbs(lbs(i)+2,i)

       b = bi * bi * z/2

       write(iout,'(a4,a4,a,3D20.8)') EBS(i),EBS(i),': ', a, b, (a-b)/a

      End do

      Close(iout)

    END SUBROUTINE test_azl




!========================================================================
      SUBROUTINE test_vc(z)
!========================================================================

! ... tests the accuracy of AZL subprogram on the basis
! ... of hydrogenic orbitals. Results are output in file 'azl_out'.

      USE spline_param
      USE spline_galerkin
      USE spline_grid
      USE spline_orbitals
      USE spline_hl

      IMPLICIT NONE

      REAL(KIND=8), INTENT(in) :: Z
      REAL(KIND=8), EXTERNAL :: HNORM, BVMV, BVALU2, AZL
      INTEGER :: iout, i,j, l,ns1, ii, jj
      REAL(KIND=8) :: a,b,ai,bi,aj,bj, b2s, b2ss, b3ss,fine
      REAL(KIND=8), DIMENSION(nv,ks,ks) :: dbiatx
      REAL(KIND=8), DIMENSION(nv) :: x,y,y1,y2
      REAL(KIND=8), DIMENSION(ns) :: bcoef, eval
      REAL(KIND=8), DIMENSION(ns,ns) :: hm,c


      INTEGER :: Lwork, info
      REAL(KIND=8), DIMENSION(3*ns) :: work
      Lwork = 3*ns

      iout=2
      Open(iout,file='vc_out')


      write(iout,'(/a/)') '   accuracy of mass-velocity terms: '

                          !   - 0.5 * (Z/n)^4 * [ 4n/(l+1/2) -3 ]
      Do i=1,nbf

       a = lbs(i)+0.5
       a = 4*nbs(i)/a - 3
       a = -0.5 * (Z/nbs(i))**4 * a

       Call mvc(lbs(i))

       b = -0.5 * BVMV (ns,ks,vc,'s',pbs(1,i),pbs(1,i))

       write(iout,'(a4,3D20.8)') EBS(i), a,b, (a-b)/a

      End do

! ... coulomb problem ...


    Do l = 0,3

     Call HLM(l)

     fine = 0.25D0/(137.036D0)**2

     Call mvc(l)
     hl = hl + fine*vc
     if(l.eq.0) then
      a = azl(z,h,ks,l+1)
      hl(2,ks) = hl(2,ks) - z*a*a*fine
     end if

      hm = 0.d0
      c = 0.d0

       Do j = 1,ks
        Do i = ks-j+1,ns
         ii=i      -1 - l
         jj=i-ks+j -1 - l
         if(ii.lt.1.or.jj.lt.1) Cycle
         hm(ii,jj) = -0.5d0*hl(i,j)
         hm(jj,ii) = -0.5d0*hl(i,j)
         c(ii,jj) =  sb(i,j)
         c(jj,ii) =  sb(i,j)
        End do
       End do

      ns1 = ns  - l - 3
      Call DSYGV(1,'V','U',ns1,hm,ns,C,ns,eval,WORK,LWORK,INFO)

      write(iout,'(a,i5)') ' l =  ', l
      write(iout,'(5E12.5)') Eval(1:ns1)

     End do


         x = 0.0
         call vbsplvd(t,ks,1,x(1),3,dbiatx)

         write(iout,'(/a,f10.5/)')  't = ', x(1)
         write(iout,'(8E12.3)') dbiatx(1,1:ks,1)
         write(iout,'(8E12.3)') dbiatx(1,1:ks,2)
         write(iout,'(8E12.3)') dbiatx(1,1:ks,3)

         b2s = dbiatx(1,2,2)
         b2ss = dbiatx(1,2,3)
         b3ss = dbiatx(1,3,3)

      write(iout,'(/a/)') ' vc - matrix : '

      Call mvc(0)

!      vc(2,ks) = vc(2,ks) + b2s*b2ss
!      vc(3,ks-1) = vc(3,ks-1) + 0.5*b2s*b3ss

      Do i = 1,ks
       VC(i,ks-i+1) = 0.d0
      End do


      Do i = 1,ks
       write(iout,'(8E10.3)') vc(i,1:ks)
      End do


      Do i = 1, ks

       bcoef(:) = 0.d0
       bcoef(i) = 1.d0

       Do j = 1, ks
        x(j) = t(ks+j-1)
        if(j.gt.1) x(j) = x(j) - 1.d-20
        y (j) = bvalu2 (t, bcoef, ns, ks, x(j), 0)
        y1(j) = bvalu2 (t, bcoef, ns, ks, x(j), 1)
        y2(j) = bvalu2 (t, bcoef, ns, ks, x(j), 2)
       End do

!       write(iout,'(a,i5)') ' spline ', i
!       write(iout,'(8E12.3)') y(1:ks)
!       write(iout,'(8E12.3)') y1(1:ks)
!       write(iout,'(8E12.3)') y2(1:ks)

      End do


      Close(iout)

    END SUBROUTINE test_vc




























!    below are some additional routines from my libraries




!======================================================================
    Function RRTC()
!======================================================================

!   give the current time in seconds

    IMPLICIT NONE

    REAL(KIND=8) :: RRTC
    CHARACTER(LEN=8) :: D
    CHARACTER(LEN=10) :: T
    INTEGER :: id,ih,im,is,ims
    REAL TIME(2), ETIME

    Call DATE_AND_TIME(date=D,time=T)

    read(D,'(6x,i2)') id
    read(T,'(3i2,1x,i3)') ih,im,is,ims
    RRTC = id*86400 + ih*3600 + im*60 + is
    RRTC = RRTC + ims/1000.d0
    RRTC = ETIME(TIME)
    RRTC = TIME(1)*1000

    End Function RRTC


!======================================================================
      Subroutine EL4_NLK(EL,n,l,k)
!======================================================================
!
!     decodes the specroscopic notation for electron orbital (n,l,k)
!
!     It is allowed following notations: 1s, 2s3, 2p30, 20p3,
!     1sh, 20sh, kp, kp1, kp11, ns, ns3, ns33, ... .
!======================================================================

      IMPLICIT NONE
!
      CHARACTER(LEN=4), INTENT(in) :: EL
      INTEGER, INTENT(out) :: n,l,k
!
      INTEGER :: i,j,ic, LA
!
      n=0
      l=-1
      k=0
      i=1
      j=1
    1 if(EL(i:i).eq.' ') then
       i=i+1
       if(i.le.4) go to 1
      end if
      if(i.eq.4.and.j.eq.1) j=2
      if(i.gt.4) Return

      ic=ICHAR(EL(i:i))-ICHAR('1')+1

      if(j.eq.1) then

       if(n.eq.0.and.ic.le.9) then
        n=ic
       elseif(n.eq.0.and.ic.gt.9) then
        n=ic
        j=2
       elseif(ic.gt.9) then
        l=LA(EL(i:i))
        j=3
       else
        n=n*10+ic
       end if

      elseif(j.eq.2) then
       l=LA(EL(i:i))
       j=3

      elseif(j.eq.3) then
       if(ic.gt.9.and.k.eq.0) then
        k=ic
        j=4
       elseif(ic.gt.9.and.k.gt.0) then
        go to 99
       else
        k=k*10+ic
       end if
      end if

      i=i+1
      if(i.le.4.and.j.le.3) go to 1
      if(n.ge.100.or.l.lt.0.or.k.ge.100) go to 99
      Return

   99 write(*,*) ' EL4_nlk is fail to decode the ',  &
                 ' specroscopic notation: ',EL
      Stop

      End Subroutine EL4_NLK



!======================================================================
      Function ELF4(n,l,k)
!======================================================================
!     gives the specroscopic notation for electron orbital (n,l,k)
!
!     n and k must be < 100; if they > 50, they are considered
!     as characters with code i+ICHAR(1)-1.
!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(in) :: n,l,k
      CHARACTER(len=4) :: ELF4, EL
      CHARACTER(len=1) :: AL
      INTEGER :: i,kk,i1
!
      EL='    '
      i=4
      i1 = ICHAR('1')
!
      if(k.gt.0) then
       kk=k+i1-1
       if(kk.ge.97.and.kk.le.122) then           ! from a to z
        EL(i:i)=CHAR(kk)
        i=i-1
       elseif(k.lt.10) then
        write(EL(i:i),'(i1)') k
        i=i-1
       elseif(k.ge.10.and.k.lt.100) then
        write(EL(i-1:i),'(i2)') k
        i=i-2
       else
        write(*,*) ' ELF4: k is out of limits:',k
        Stop
       end if
      end if
!
      EL(i:i)=AL(l,1)
      i=i-1
!
      if(n.gt.0) then
       kk=n+i1-1
       if(kk.ge.97.and.kk.le.122) then           ! from a to z
        EL(i:i)=CHAR(kk)
       elseif(n.lt.10) then
        write(EL(i:i),'(i1)') n
       elseif(n.ge.10.and.n.lt.100.and.i.ge.2) then
        write(EL(i-1:i),'(i2)') n
       else
        write(*,*) ' ELF4: n is out of limits:',n
        Stop
       end if
      end if
!
      ELF4=EL

      END FUNCTION ELF4



!====================================================================
      Function LA(a)
!====================================================================

      Character a, SET*20

      Data SET/'spdfghiklmSPDFGHIKLM'/

      i=1000
      DO j=1,20
       if(a.eq.SET(j:j)) then
        i=j
        Exit
       end if
      End do

      if(i.le.10) then
       la=i-1
      elseif(i.le.20) then
       la=i-11
      else
       la=10
      end if
       if(la.lt.0) la=10

      End function LA


!====================================================================
      CHARACTER FUNCTION AL(L,K)
!====================================================================

! ..  provides some spectroscopic symbols

      CHARACTER AS*11,AB*11,AN*16
      DATA AS/'spdfghiklm*'/,       &
           AB/'SPDFGHIKLM*'/,       &
           AN/'0123456789ABCDEF'/

      I=L+1
      IF(K.EQ.5.OR.K.EQ.6) I=(L-1)/2+1
      AL='?'
      IF(I.GE.1.AND.I.LE.11) THEN
        IF(K.EQ.1) AL=AS(I:I)
        IF(K.EQ.2) AL=AB(I:I)
        IF(K.EQ.3) AL=AN(I:I)
        IF(K.EQ.5) AL=AS(I:I)
        IF(K.EQ.6) AL=AB(I:I)
      END IF
      IF(K.EQ.4) THEN
        if(L.eq.-1.OR.L.eq.0) AL='o'
        if(L.eq.+1          ) AL='e'
      END IF
      IF(K.EQ.7) THEN
       if(L.le.0) AL='-'
       if(L.gt.0) AL='+'
       END IF
      END FUNCTION AL






!======================================================================
      FUNCTION rkv(i1,j1,i2,j2,k)
!======================================================================
!
!                 k
!     Evaluates  N (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    REAL(KIND=8), DIMENSION(ns) :: v
    REAL(KIND=8) :: rkv
    INTEGER :: icase, int

      rkv = 0.d0

      Call Mrk_cell(k)

      Call INT_v(i1,j1,i2,j2,k,5,1,v)
      rkv =  rkv + SUM(v(:)*p(:,i1))

      Call INT_v(i1,j1,i2,j2,k,5,2,v)
      rkv =  rkv + SUM(v(:)*p(:,j1))

      Call INT_v(i1,j1,i2,j2,k,5,3,v)
      rkv =  rkv + SUM(v(:)*p(:,i2))

      Call INT_v(i1,j1,i2,j2,k,5,4,v)
      rkv =  rkv + SUM(v(:)*p(:,j2))

      rkv = rkv * fine / 4

      END FUNCTION rkv



!======================================================================
    FUNCTION nkd(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  N (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: a,b
    REAL(KIND=8) :: nkd
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    nkd = 0.d0

    Call Mnk_cell(k)

    Call INT_de (p(1,i1),p(1,i2),a,8,2,sym)
    Call density(ns,ks,b,p(1,j1),p(1,j2),sym)

    nkd = nkd + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,j1),p(1,j2),a,8,1,sym)
    Call density(ns,ks,b,p(1,i1),p(1,i2),sym)

    nkd = nkd + SUM_amb(ns,ks,a,b,sym)

    nkd = nkd * fine /2

    END FUNCTION nkd



!======================================================================
    FUNCTION nke(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  N (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    INTEGER :: i,j
    REAL(KIND=8), DIMENSION(ns,ns) :: a,b
    REAL(KIND=8) :: nke
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    nke = 0.d0

    Call Mnk_cell(k)

    Call INT_de (p(1,i2),p(1,j1),a,8,3,sym)
    Call density(ns,ks,b,p(1,i1),p(1,j2),sym)
    nke = nke + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,i1),p(1,j2),a,8,4,sym)
    Call density(ns,ks,b,p(1,i2),p(1,j1),sym)
    nke = nke + SUM_amb(ns,ks,a,b,sym)

    nke = nke * fine / 2

    END FUNCTION nke


!======================================================================
      FUNCTION nkv(i1,j1,i2,j2,k)
!======================================================================
!
!                 k
!     Evaluates  N (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    REAL(KIND=8), DIMENSION(ns) :: v
    REAL(KIND=8) :: nkv
    INTEGER :: icase, int

      nkv = 0.d0

      Call Mnk_cell(k)

      Call INT_v(i1,j1,i2,j2,k,8,1,v)
      nkv =  nkv + SUM(v(:)*p(:,i1))

      Call INT_v(i1,j1,i2,j2,k,8,2,v)
      nkv =  nkv + SUM(v(:)*p(:,j1))

      Call INT_v(i1,j1,i2,j2,k,8,3,v)
      nkv =  nkv + SUM(v(:)*p(:,i2))

      Call INT_v(i1,j1,i2,j2,k,8,4,v)
      nkv =  nkv + SUM(v(:)*p(:,j2))

      nkv = nkv * fine / 4

      END FUNCTION nkv




!======================================================================
    FUNCTION vkd(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  V (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: a,b
    REAL(KIND=8) :: vkd
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    vkd = 0.d0

    Call Mvk_cell(k)

    Call INT_de (p(1,i1),p(1,i2),a,9,2,sym)
    Call density(ns,ks,b,p(1,j1),p(1,j2),sym)

    vkd = vkd + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,j1),p(1,j2),a,9,1,sym)
    Call density(ns,ks,b,p(1,i1),p(1,i2),sym)

    vkd = vkd + SUM_amb(ns,ks,a,b,sym)

    vkd = vkd * fine /2

    END FUNCTION vkd



!======================================================================
    FUNCTION vke(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  V (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    INTEGER :: i,j
    REAL(KIND=8), DIMENSION(ns,ns) :: a,b
    REAL(KIND=8) :: vke
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    vke = 0.d0

    Call Mvk_cell(k)

    Call INT_de (p(1,i2),p(1,j1),a,9,3,sym)
    Call density(ns,ks,b,p(1,i1),p(1,j2),sym)
    vke = vke + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,i1),p(1,j2),a,9,4,sym)
    Call density(ns,ks,b,p(1,i2),p(1,j1),sym)
    vke = vke + SUM_amb(ns,ks,a,b,sym)

    vke = vke * fine / 2

    END FUNCTION vke


!======================================================================
      FUNCTION vkv(i1,j1,i2,j2,k)
!======================================================================
!
!                 k
!     Evaluates  V (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    REAL(KIND=8), DIMENSION(ns) :: v
    REAL(KIND=8) :: vkv
    INTEGER :: icase, int

      vkv = 0.d0

      Call Mvk_cell(k)

      Call INT_v(i1,j1,i2,j2,k,9,1,v)
      vkv =  vkv + SUM(v(:)*p(:,i1))

      Call INT_v(i1,j1,i2,j2,k,9,2,v)
      vkv =  vkv + SUM(v(:)*p(:,j1))

      Call INT_v(i1,j1,i2,j2,k,9,3,v)
      vkv =  vkv + SUM(v(:)*p(:,i2))

      Call INT_v(i1,j1,i2,j2,k,9,4,v)
      vkv =  vkv + SUM(v(:)*p(:,j2))

      vkv = vkv * fine / 4

      END FUNCTION vkv






!======================================================================
    FUNCTION mkd(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  M (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: a,b
    REAL(KIND=8) :: mkd
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    mkd = 0.d0

    Call Mmk_cell(k)

    Call INT_de (p(1,i1),p(1,i2),a,4,2,sym)
    Call density(ns,ks,b,p(1,j1),p(1,j2),sym)

    mkd = mkd + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,j1),p(1,j2),a,4,1,sym)
    Call density(ns,ks,b,p(1,i1),p(1,i2),sym)

    mkd = mkd + SUM_amb(ns,ks,a,b,sym)

    mkd = mkd * fine /2  /2

    END FUNCTION mkd



!======================================================================
    FUNCTION mke(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  M (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    INTEGER :: i,j
    REAL(KIND=8), DIMENSION(ns,ns) :: a,b
    REAL(KIND=8) :: mke
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    mke = 0.d0

    Call Mmk_cell(k)

    Call INT_de (p(1,i2),p(1,j1),a,4,3,sym)
    Call density(ns,ks,b,p(1,i1),p(1,j2),sym)
    mke = mke + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,i1),p(1,j2),a,4,4,sym)
    Call density(ns,ks,b,p(1,i2),p(1,j1),sym)
    mke = mke + SUM_amb(ns,ks,a,b,sym)

    mke = mke * fine / 2  /2

    END FUNCTION mke


!======================================================================
      FUNCTION mkv(i1,j1,i2,j2,k)
!======================================================================
!
!                 k
!     Evaluates  M (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    REAL(KIND=8), DIMENSION(ns) :: v
    REAL(KIND=8) :: mkv
    INTEGER :: icase, int

      mkv = 0.d0

      Call Mmk_cell(k)

      Call INT_v(i1,j1,i2,j2,k,4,1,v)
      mkv =  mkv + SUM(v(:)*p(:,i1))

      Call INT_v(i1,j1,i2,j2,k,4,2,v)
      mkv =  mkv + SUM(v(:)*p(:,j1))

      Call INT_v(i1,j1,i2,j2,k,4,3,v)
      mkv =  mkv + SUM(v(:)*p(:,i2))

      Call INT_v(i1,j1,i2,j2,k,4,4,v)
      mkv =  mkv + SUM(v(:)*p(:,j2))

      mkv = mkv * fine / 4  /2

      END FUNCTION mkv


!======================================================================
    FUNCTION tkd(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  T (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: a,b
    REAL(KIND=8) :: tkd
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    tkd = 0.d0

    Call Mtk_cell(k)

    Call INT_de (p(1,i1),p(1,i2),a,3,2,sym)
    Call density(ns,ks,b,p(1,j1),p(1,j2),sym)

    tkd = tkd + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,j1),p(1,j2),a,3,1,sym)
    Call density(ns,ks,b,p(1,i1),p(1,i2),sym)

    tkd = tkd + SUM_amb(ns,ks,a,b,sym)

    tkd = tkd  / 2

    END FUNCTION tkd



!======================================================================
    FUNCTION tke(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  T (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    INTEGER :: i,j
    REAL(KIND=8), DIMENSION(ns,ns) :: a,b
    REAL(KIND=8) :: tke
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    tke = 0.d0

    Call Mtk_cell(k)

    Call INT_de (p(1,i2),p(1,j1),a,3,3,sym)
    Call density(ns,ks,b,p(1,i1),p(1,j2),sym)
    tke = tke + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,i1),p(1,j2),a,3,4,sym)
    Call density(ns,ks,b,p(1,i2),p(1,j1),sym)
    tke = tke + SUM_amb(ns,ks,a,b,sym)

    tke = tke  / 2

    END FUNCTION tke


!======================================================================
      FUNCTION tkv(i1,j1,i2,j2,k)
!======================================================================
!
!                 k
!     Evaluates  T (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    REAL(KIND=8), DIMENSION(ns) :: v
    REAL(KIND=8) :: tkv
    INTEGER :: icase, int

      tkv = 0.d0

      Call Mtk_cell(k)

      Call INT_v(i1,j1,i2,j2,k,3,1,v)
      tkv =  tkv + SUM(v(:)*p(:,i1))

      Call INT_v(i1,j1,i2,j2,k,3,2,v)
      tkv =  tkv + SUM(v(:)*p(:,j1))

      Call INT_v(i1,j1,i2,j2,k,3,3,v)
      tkv =  tkv + SUM(v(:)*p(:,i2))

      Call INT_v(i1,j1,i2,j2,k,3,4,v)
      tkv =  tkv + SUM(v(:)*p(:,j2))

      tkv = tkv  / 4

      END FUNCTION tkv




!======================================================================
    FUNCTION qkd(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  Q (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: a,b
    REAL(KIND=8) :: qkd
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    qkd = 0.d0

    Call Mqk_cell(k)

    Call INT_de (p(1,i1),p(1,i2),a,2,2,sym)
    Call density(ns,ks,b,p(1,j1),p(1,j2),sym)

    qkd = qkd + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,j1),p(1,j2),a,2,1,sym)
    Call density(ns,ks,b,p(1,i1),p(1,i2),sym)

    qkd = qkd + SUM_amb(ns,ks,a,b,sym)

    qkd = qkd / 2

    END FUNCTION qkd



!======================================================================
    FUNCTION qke(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  Q (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    INTEGER :: i,j
    REAL(KIND=8), DIMENSION(ns,ns) :: a,b
    REAL(KIND=8) :: qke
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    qke = 0.d0

    Call Mqk_cell(k)

    Call INT_de (p(1,i2),p(1,j1),a,2,3,sym)
    Call density(ns,ks,b,p(1,i1),p(1,j2),sym)
    qke = qke + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,i1),p(1,j2),a,2,4,sym)
    Call density(ns,ks,b,p(1,i2),p(1,j1),sym)
    qke = qke + SUM_amb(ns,ks,a,b,sym)

    qke = qke / 2

    END FUNCTION qke


!======================================================================
      FUNCTION qkv(i1,j1,i2,j2,k)
!======================================================================
!
!                 k
!     Evaluates  Q (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    REAL(KIND=8), DIMENSION(ns) :: v
    REAL(KIND=8) :: qkv
    INTEGER :: icase, int

      qkv = 0.d0

      Call Mqk_cell(k)

      Call INT_v(i1,j1,i2,j2,k,2,1,v)
      qkv =  qkv + SUM(v(:)*p(:,i1))

      Call INT_v(i1,j1,i2,j2,k,2,2,v)
      qkv =  qkv + SUM(v(:)*p(:,j1))

      Call INT_v(i1,j1,i2,j2,k,2,3,v)
      qkv =  qkv + SUM(v(:)*p(:,i2))

      Call INT_v(i1,j1,i2,j2,k,2,4,v)
      qkv =  qkv + SUM(v(:)*p(:,j2))

      qkv = qkv / 4

      END FUNCTION qkv



!======================================================================
    FUNCTION wkd(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  W (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: a,b
    REAL(KIND=8) :: wkd
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    wkd = 0.d0

    Call Mwk_cell(k)

    Call INT_de (p(1,i1),p(1,i2),a,9,2,sym)
    Call density(ns,ks,b,p(1,j1),p(1,j2),sym)

    wkd = wkd + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,j1),p(1,j2),a,9,1,sym)
    Call density(ns,ks,b,p(1,i1),p(1,i2),sym)

    wkd = wkd + SUM_amb(ns,ks,a,b,sym)

    wkd = wkd * fine /2

    END FUNCTION wkd



!======================================================================
    FUNCTION wke(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  W (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    INTEGER :: i,j
    REAL(KIND=8), DIMENSION(ns,ns) :: a,b
    REAL(KIND=8) :: wke
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    wke = 0.d0

    Call Mwk_cell(k)

    Call INT_de (p(1,i2),p(1,j1),a,9,3,sym)
    Call density(ns,ks,b,p(1,i1),p(1,j2),sym)
    wke = wke + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,i1),p(1,j2),a,9,4,sym)
    Call density(ns,ks,b,p(1,i2),p(1,j1),sym)
    wke = wke + SUM_amb(ns,ks,a,b,sym)

    wke = wke * fine / 2

    END FUNCTION wke


!======================================================================
      FUNCTION wkv(i1,j1,i2,j2,k)
!======================================================================
!
!                 k
!     Evaluates  W (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    REAL(KIND=8), DIMENSION(ns) :: v
    REAL(KIND=8) :: wkv
    INTEGER :: icase, int

      wkv = 0.d0

      Call Mwk_cell(k)

      Call INT_v(i1,j1,i2,j2,k,9,1,v)
      wkv =  wkv + SUM(v(:)*p(:,i1))

      Call INT_v(i1,j1,i2,j2,k,9,2,v)
      wkv =  wkv + SUM(v(:)*p(:,j1))

      Call INT_v(i1,j1,i2,j2,k,9,3,v)
      wkv =  wkv + SUM(v(:)*p(:,i2))

      Call INT_v(i1,j1,i2,j2,k,9,4,v)
      wkv =  wkv + SUM(v(:)*p(:,j2))

      wkv = wkv * fine / 4

      END FUNCTION wkv


!====================================================================
    SUBROUTINE mvc(l)
!====================================================================
!
!   Computes the matrix elements for the mass-velocity correction
!   in the B-spline basis.  The lack of symmetry in the d^2/dr^2
!   operator is ignored.
!
!     VC(i,j) = INT [  (d^2/dr^2 - l(l+1)/r) B_i(r) *
!                      (d^2/dr^2 - l(l+1)/r) B_j(r)  ] dr
!--------------------------------------------------------------------
!
!   on entry
!   --------
!       l    the angular momentum
!
!   on exit
!   -------
!       vc   the mass velocity correction in symmetric storage mode
!
!--------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_atomic
    USE spline_hl

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l

    ! .. local variables

    INTEGER :: m, ith, jth, i, irow, jcol
    REAL(KIND=8) :: fll, y1, y2, S, B

    ! .. initialize the vc array

    vc = 0.d0
    fll  =  l*(l+1)

    ! .. compute the matrix elements

    do m = 1,ks
      do i = 1,nv
        S = fll*grm(i,m)*grm(i,m)

! ... cutoff correction

        B = gr(i,m)/(gr(i,m)+2*fine*Z)
        B = B*B*B

        do ith = 1,ks
          irow = i+ith-1
          do jth = 1,ith
          jcol = jth-ith+ks

            y1 = bspd(i,m,ith,2) - S*bsp(i,m,ith)
            y2 = bspd(i,m,jth,2) - S*bsp(i,m,jth)
            vc(irow,jcol) = vc(irow,jcol) + grw(i,m)*y1*y2 *B

          end do
        end do
      end do
    end do

    END SUBROUTINE mvc



!======================================================================
      REAL(8) FUNCTION vkc(i1,j1,i2,j2,k)
!======================================================================
!               k
!     Returns  V (i1, j1; i2, j2) - direct summations of moments
!                                   over cells
!----------------------------------------------------------------------

      USE spline_param
      USE spline_orbitals, p => pbs
      USE spline_moments
      USE spline_atomic

      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: i1,j1,i2,j2,k

      ! .. local variables

      INTEGER(4) :: i,j, iv, ivi, ivj, ii, ik,jk
      REAL(8), DIMENSION(nv) :: v1, v2, v3, v4
      REAL(8), DIMENSION(ks*ks) ::  a, b
      REAL(8) :: s1, s2

      ! .. check the need of calculations

      Call vk_moments(k)

      ik = ks*(ks+1)/2
      jk = ks*ks

      vkc = 0.d0

      Do iv = 1,nv

       ii=0
       Do i=1,ks
        ivi=iv+i-1
        Do j=1,ks
         ivj=iv+j-1
         ii = ii+1
         a(ii) = p(ivi,i1)*p(ivj,i2)
        End do
       End do

       ii=0
       Do i=1,ks
        ivi=iv+i-1
        Do j=i,ks
         ivj=iv+j-1
         ii = ii+1
         if(i.eq.j) then
          b(ii) = p(ivi,j1)*p(ivj,j2)
         else
          b(ii) = p(ivi,j1)*p(ivj,j2) + p(ivi,j2)*p(ivj,j1)
         end if
        End do
       End do

       v1(iv) = SUM(rkd3(1:jk,iv)*a)
       v2(iv) = SUM(rkd2(1:ik,iv)*b)
       v3(iv) = SUM(rkd4(1:jk,iv)*a)
       v4(iv) = SUM(rkd1(1:ik,iv)*b)

       ! the diagonal cell contribution

       Do j = 1,ik
        vkc = vkc + SUM(a(1:jk)*rkd(1:jk,j,iv))*b(j)
       End do

      End do

      ! the upper and lower regions

      s1 = 0.d0
      s2 = 0.d0
      Do iv =  2,nv
       s1 = s1 + v1(iv-1)
       vkc = vkc + s1*v2(iv)
       s2 = s2 + v4(iv-1)
       vkc = vkc + s2*v3(iv)
      End do

      vkc = vkc * fine

      End FUNCTION vkc


!======================================================================
    REAL(8) FUNCTION rkdn(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  r (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_densities
    USE spline_Rk_integrals

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    REAL(8), External :: SUM_amb

    Call make_rkd(k,i1,i2)
    Call make_density(j1,j2,'s')

    rkdn = SUM_amb(ns,ks,ds,rkbd,'l')

    Call make_rkd(k,j1,j2)
    Call make_density(i1,i2,'s')

    rkdn = rkdn + SUM_amb(ns,ks,ds,rkbd,'l')

    rkdn = rkdn / 2.d0

    END FUNCTION rkdn

!======================================================================
    REAL(8) FUNCTION rken(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  r (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_densities
    USE spline_Rk_integrals

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    REAL(8), External :: SUM_amb

    Call make_rke(k,j1,i2)
    Call make_density(i1,j2,'x')

    rken = SUM_amb(ns,ks,dx,rkbe,'x')

    Call make_rke(k,i1,j2)
    Call make_density(j1,i2,'x')

    rken = rken + SUM_amb(ns,ks,dx,rkbe,'x')

    rken = rken / 2.d0

    END FUNCTION rken

!======================================================================
    REAL(8) FUNCTION rkvn(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  r (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_densities
    USE spline_Rk_integrals
    Use spline_orbitals

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    Call make_rkbv(k,i1,j1,i2)
    rkvn = SUM(pbs(:,j2)*rkbv(:))

    Call make_rkbv(k,i1,j2,i2)
    rkvn = rkvn + SUM(pbs(:,j1)*rkbv(:))

    Call make_rkbv(k,j1,i1,j2)
    rkvn = rkvn + SUM(pbs(:,i2)*rkbv(:))

    Call make_rkbv(k,j1,i2,j2)
    rkvn = rkvn + SUM(pbs(:,i1)*rkbv(:))

    rkvn = rkvn / 4.d0

    END FUNCTION rkvn


!======================================================================
    REAL(8) FUNCTION rk1n(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  r (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_densities
    USE spline_Rk_integrals
    Use spline_orbitals

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k
    REAL(8) :: d(ns,ks), v(ns)
    Integer(4) :: i
    REAL(8), External :: SUM_amb

    Call make_rkb1(k,i1)
    Call make_density(j1,j2,'s')
    Do i=1,ns
     d(1:ns,1:ks) = rkb1(i,1:ns,1:ks)
     v(i) = SUM_amb(ns,ks,ds,d,'l')
    End do
    rk1n = SUM(pbs(:,i2)*v(:))

    Call make_rkb1(k,i2)
    Call make_density(j1,j2,'s')
    Do i=1,ns
     d(1:ns,1:ks) = rkb1(i,1:ns,1:ks)
     v(i) = SUM_amb(ns,ks,ds,d,'l')
    End do
    rk1n = rk1n + SUM(pbs(:,i1)*v(:))

    Call make_rkb1(k,j1)
    Call make_density(i1,i2,'s')
    Do i=1,ns
     d(1:ns,1:ks) = rkb1(i,1:ns,1:ks)
     v(i) = SUM_amb(ns,ks,ds,d,'l')
    End do
    rk1n = rk1n + SUM(pbs(:,j2)*v(:))

    Call make_rkb1(k,j2)
    Call make_density(i1,i2,'s')
    Do i=1,ns
     d(1:ns,1:ks) = rkb1(i,1:ns,1:ks)
     v(i) = SUM_amb(ns,ks,ds,d,'l')
    End do
    rk1n = rk1n + SUM(pbs(:,j1)*v(:))

    rk1n = rk1n / 4.d0

    END FUNCTION rk1n


!======================================================================
    FUNCTION rkd(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  N (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: a,b
    REAL(KIND=8) :: rkd
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    rkd = 0.d0

    Call Mrk_cell(k)

    Call INT_de (p(1,i1),p(1,i2),a,5,2,sym)
    Call density(ns,ks,b,p(1,j1),p(1,j2),sym)

    rkd = rkd + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,j1),p(1,j2),a,5,1,sym)
    Call density(ns,ks,b,p(1,i1),p(1,i2),sym)

    rkd = rkd + SUM_amb(ns,ks,a,b,sym)

    rkd = rkd * fine /2

    END FUNCTION rkd


!======================================================================
    FUNCTION rke(i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  N (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    INTEGER :: i,j
    REAL(KIND=8), DIMENSION(ns,ns) :: a,b
    REAL(KIND=8) :: rke
    Character(1) :: sym
    REAL(8), External :: SUM_amb

    rke = 0.d0

    Call Mrk_cell(k)

    Call INT_de (p(1,i2),p(1,j1),a,5,3,sym)
    Call density(ns,ks,b,p(1,i1),p(1,j2),sym)
    rke = rke + SUM_amb(ns,ks,a,b,sym)

    Call INT_de (p(1,i1),p(1,j2),a,5,4,sym)
    Call density(ns,ks,b,p(1,i2),p(1,j1),sym)
    rke = rke + SUM_amb(ns,ks,a,b,sym)

    rke = rke / 2

    END FUNCTION rke
