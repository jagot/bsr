!======================================================================
      PROGRAM test_int
!======================================================================
!     Testing two-electron integrals: list of desirable integrals are
!     supposed to be in name.int file, the corresponding one-electron
!     orbitals - in name.bsw.
!     Results - in name.test
!     B-splines requires knot.dat file
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_nuclear, z => atomic_number
      Use DBS_gauss
      Use DBS_orbitals_pq

      Implicit real(8) (A-H,O-Z)

      Character(80) :: name,AF
      Integer :: nui=2; Character(80) :: AF_int = 'int.inp'
      Integer :: nut=3; Character(80) :: AF_out = 'int.out'

! ... define B-splines:


      zz = 0.d0; aa = 2.0
      Call def_grid('knot.dat',' ',zz,aa)

      Call alloc_DBS_gauss

! ... open files:

      Call Read_aarg('AF_int',AF_int)
      Call Read_aarg('AF_out',AF_out)

! ... generate

      nn=3;  Call Read_iarg('nn',nn)

      Do n=1,nn
       Do k = -n,n-1
        if(k.eq.0) Cycle
        i = Ifind_bsorb(n,k,0,2)
        CALL bdcwf_pq(n,k,z,pq(1,1,i),pq(1,2,i)); mbs(i)=ns
       End do
      End do
      kmax = maxval(jbs(1:nbf))+1

      Call alloc_RK_integrals(ns,ks,0,0,4)
!      Call Print_cwf_pq

      Open(nui,file=AF_int)
      Open(nut,file=AF_out)

      write(nut,'(a,i5)') 'ksp =',ksp
      write(nut,'(a,i5)') 'ksq =',ksq
      write(nut,*)

      write(nut,*) 'number of radial functions:',nbf
      write(nut,'(10a5)') EBS(1:nbf)
      write(nut,*)

      CALL test_integrals(nui,nut)

      End Program test_int



!========================================================================
      SUBROUTINE test_integrals(nui,nut)
!========================================================================
!     tests the given set of radial integrals
!     The integrals are given in unit nui, results - in unit nut.
!------------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)

      Character(1) :: INT
      Character(80) :: AS
      Character( 5) :: EL1, EL2, EL3, EL4

      Integer :: i,k, in,iout

      Integer :: n1,n2,n3,n4
      Integer :: l1,l2,l3,l4
      Integer :: k1,k2,k3,k4
      Integer :: i1,i2,i3,i4
      Integer :: ie1,ie2,ie3,ie4

      Real(8) :: S, SS, SA,SB, S1, S2

      Real(8) :: t1,t2

      Integer, external :: Ifind_bsorb
      Real(8), external :: RRTC

!----------------------------------------------------------------------
! ... loop over input strings

    1 read(nui,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2          ! end of data
      if(LEN_TRIM(AS).eq.0) go to 1
!----------------------------------------------------------------------
      if(AS(1:1).eq.'#') then             ! timing block
       go to 1
      end if                              ! end of timing mode (#)

!----------------------------------------------------------------------
! ... test of accuracy:

      Read(AS,'(a1,i2,4(1x,a5))') int,k,el1,el2,el3,el4
      m=0
      if(int.eq.'R') m=1
      if(int.eq.'S') m=2
      if(m.eq.0) go to 1

      Call EL_NLJK(el1,n1,k1,l1,j1,i1);  ie1=Ifind_bsorb(n1,k1,i1,0)
      if(ie1.eq.0) go to 1
      Call EL_NLJK(el2,n2,k2,l2,j2,i2);  ie2=Ifind_bsorb(n2,k2,i2,0)
      if(ie2.eq.0) go to 1
      Call EL_NLJK(el3,n3,k3,l3,j3,i3);  ie3=Ifind_bsorb(n3,k3,i3,0)
      if(ie3.eq.0) go to 1
      Call EL_NLJK(el4,n4,k4,l4,j4,i4);  ie4=Ifind_bsorb(n4,k4,i4,0)
      if(ie4.eq.0) go to 1

      i = INDEX(AS,'='); if(i.eq.0) go to 1
      read(AS(i+1:),*) SA

      write(nut,'(a,d24.16)') AS(1:30),SA
      write(*  ,'(a,d24.16)') AS(1:30),SA

!----------------------------------------------------------------------
      if(int.eq.'R') then

       S = rk1 (ie1,ie2,ie3,ie4,k);  if(SA.eq.0) SA=S
       SS = (S-SA)/SA
       write(nut,'(30x,d24.16,d12.3,5x,a)') s,ss,'  rk1(.,a;.,b)'
       write(*  ,'(30x,d24.16,d12.3,5x,a)') s,ss,'  rk1(.,a;.,b)'
       S = rk2 (ie1,ie2,ie3,ie4,k)
       SS = (S-SA)/SA
       write(nut,'(30x,d24.16,d12.3,5x,a)') s,ss,'  rk2(a,.;b,.)'
       write(*  ,'(30x,d24.16,d12.3,5x,a)') s,ss,'  rk2(a,.;b,.)'
       S = rk3 (ie1,ie2,ie3,ie4,k)
       SS = (S-SA)/SA
       write(nut,'(30x,d24.16,d12.3,5x,a)') s,ss,'  rk3(.,a;b,.)'
       write(*  ,'(30x,d24.16,d12.3,5x,a)') s,ss,'  rk3(.,a;b,.)'
       S = rk4 (ie1,ie2,ie3,ie4,k)
       SS = (S-SA)/SA
       write(nut,'(30x,d24.16,d12.3,5x,a)') s,ss,'  rk4(a,.;.,b)'
       write(*  ,'(30x,d24.16,d12.3,5x,a)') s,ss,'  rk4(a,.;..b)'
       S = hf_gk1 (ie1,ie2,k)
       SS = (S-SA)/SA
       write(nut,'(30x,d24.16,d12.3,5x,a)') s,ss,'  gk4(a,.;.,b)'
       write(*  ,'(30x,d24.16,d12.3,5x,a)') s,ss,'  gk4(a,.;..b)'

      end if

!----------------------------------------------------------------------

      if(INT.eq.'S') then

       S = sk_ppqq (ie1,ie2,ie3,ie4,k); if(sa.eq.0.d0) SA=s

       SS = (S-SA)/SA

       write(nut,'(30x,d24.16,d12.3,5x,a)') s,ss,'  sk_ppqq'
       write(*  ,'(30x,d24.16,d12.3,5x,a)') s,ss,'  sk_ppqq'

      end if

      go to 1    ! go to new integral

    2 Continue   ! end of data

      End Subroutine test_integrals


!======================================================================
      Real(8) Function rk1 (i1,j1,i2,j2,k)
!======================================================================
!     INT = 1  - convolution other 2 and 4 variables, RK(.a;.b)
!----------------------------------------------------------------------
      Use DBS_grid,         only: ns,ks
      Use DBS_orbitals_pq,  only: p => pq

      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens1(ns,ks),dens2(ns,ks),dens3(ns,ks),dens4(ns,ks), &
                 conv(ns,ks), t1,t2
      Real(8), external :: SUM_AmB, RRTC
      Integer :: i

      rk1 = 0.d0

      Call density (ns,ks,dens1,p(1,1,i1),p(1,1,i2),'s')
      Call density (ns,ks,dens2,p(1,1,j1),p(1,1,j2),'s')
      Call density (ns,ks,dens3,p(1,2,i1),p(1,2,i2),'s')
      Call density (ns,ks,dens4,p(1,2,j1),p(1,2,j2),'s')

      Call mrk_pppp(k)
      Call convol  (ns,ks,conv,dens2,1,'s','s')
      rk1 = rk1 + SUM_AmB(ns,ks,conv,dens1,'s')
      Call mrk_qqqq(k)
      Call convol  (ns,ks,conv,dens4,1,'s','s')
      rk1 = rk1 + SUM_AmB(ns,ks,conv,dens3,'s')

      Call mrk_qpqp(k)
      Call convol  (ns,ks,conv,dens2,1,'s','s')
      rk1 = rk1 + SUM_AmB(ns,ks,conv,dens3,'s')

      Call mrk_pqpq(k)
      Call convol  (ns,ks,conv,dens4,1,'s','s')
      rk1 = rk1 + SUM_AmB(ns,ks,conv,dens1,'s')


 t1 = RRTC()
 Do i=1,10000
  Call convol  (ns,ks,conv,dens2,1,'s','s')
 End do

 t2 = RRTC()

 write(*,*) '1, s,s, 10000', t2-t1


      End Function rk1







!======================================================================
      Real(8) Function rk2 (i1,j1,i2,j2,k)
!======================================================================
!     INT = 2  - convolution other 1 and 3 variables, RK(a.;b.)
!----------------------------------------------------------------------
      Use DBS_grid,         only: ns,ks
      Use DBS_orbitals_pq,  only: p => pq

      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens1(ns,ks),dens2(ns,ks),dens3(ns,ks),dens4(ns,ks), &
                 conv(ns,ks)
      Real(8), external :: SUM_AmB

      rk2 = 0.d0

      Call density (ns,ks,dens1,p(1,1,i1),p(1,1,i2),'s')
      Call density (ns,ks,dens2,p(1,1,j1),p(1,1,j2),'s')
      Call density (ns,ks,dens3,p(1,2,i1),p(1,2,i2),'s')
      Call density (ns,ks,dens4,p(1,2,j1),p(1,2,j2),'s')

      Call mrk_pppp(k)
      Call convol  (ns,ks,conv,dens1,2,'s','s')
      rk2 = rk2 + SUM_AmB(ns,ks,conv,dens2,'s')

      Call mrk_qqqq(k)
      Call convol  (ns,ks,conv,dens3,2,'s','s')
      rk2 = rk2 + SUM_AmB(ns,ks,conv,dens4,'s')

      Call mrk_qpqp(k)
      Call convol  (ns,ks,conv,dens3,2,'s','s')
      rk2 = rk2 + SUM_AmB(ns,ks,conv,dens2,'s')

      Call mrk_pqpq(k)
      Call convol  (ns,ks,conv,dens1,2,'s','s')
      rk2 = rk2 + SUM_AmB(ns,ks,conv,dens4,'s')

      End Function rk2


!======================================================================
      Real(8) Function rk3 (i1,j1,i2,j2,k)
!======================================================================
!     INT = 3  - convolution other 2 and 3 variables, RK(.a;b.)
!----------------------------------------------------------------------
      Use DBS_grid,         only: ns,ks
      Use DBS_orbitals_pq,  only: p => pq

      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens(ns,ns), conv(ns,ns), RRTC, t1,t2
      Real(8), external :: SUM_AmB
      Integer :: i

      rk3 = 0.d0

      Call mrk_pppp(k)
      Call density (ns,ks,dens,p(1,1,i2),p(1,1,j1),'x')
      Call convol  (ns,ks,conv,dens,3,'s','s')
      Call density (ns,ks,dens,p(1,1,i1),p(1,1,j2),'x')
      rk3 = rk3 + SUM_AmB(ns,ks,conv,dens,'x')

 t1 = RRTC()
 Do i=1,100
  Call convol  (ns,ks,conv,dens,3,'s','s')
 End do

 t2 = RRTC()

 write(*,*) '3, x,x, 1000', t2-t1


      Call mrk_qqqq(k)
      Call density (ns,ks,dens,p(1,2,i2),p(1,2,j1),'x')
      Call convol  (ns,ks,conv,dens,3,'s','s')
      Call density (ns,ks,dens,p(1,2,i1),p(1,2,j2),'x')
      rk3 = rk3 + SUM_AmB(ns,ks,conv,dens,'x')

      Call mrk_qpqp(k)
      Call density (ns,ks,dens,p(1,2,i2),p(1,1,j1),'x')
      Call convol  (ns,ks,conv,dens,3,'s','s')
      Call density (ns,ks,dens,p(1,2,i1),p(1,1,j2),'x')
      rk3 = rk3 + SUM_AmB(ns,ks,conv,dens,'x')

      Call mrk_pqpq(k)
      Call density (ns,ks,dens,p(1,1,i2),p(1,2,j1),'x')
      Call convol  (ns,ks,conv,dens,3,'s','s')
      Call density (ns,ks,dens,p(1,1,i1),p(1,2,j2),'x')
      rk3 = rk3 + SUM_AmB(ns,ks,conv,dens,'x')

      End Function rk3


!======================================================================
      Real(8) Function rk4 (i1,j1,i2,j2,k)
!======================================================================
!     INT = 4  - convolution other 1 and 4 variables, RK(a.;.b)
!----------------------------------------------------------------------
      Use DBS_grid,         only: ns,ks
      Use DBS_orbitals_pq,  only: p => pq

      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens(ns,ns), conv(ns,ns),  RRTC, t1,t2
      Real(8), external :: SUM_AmB
      Integer :: i

      rk4 = 0.d0

      Call mrk_pppp(k)
      Call density (ns,ks,dens,p(1,1,i1),p(1,1,j2),'x')
      Call convol  (ns,ks,conv,dens,4,'s','s')
      Call density (ns,ks,dens,p(1,1,i2),p(1,1,j1),'x')
      rk4 = rk4 + SUM_AmB(ns,ks,conv,dens,'x')


 t1 = RRTC()
 Do i=1,1000
  Call convol  (ns,ks,conv,dens,4,'s','s')
 End do
 t2 = RRTC()

 write(*,*) '4, x,x, 1000', t2-t1





      Call mrk_qqqq(k)
      Call density (ns,ks,dens,p(1,2,i1),p(1,2,j2),'x')
      Call convol  (ns,ks,conv,dens,4,'s','s')
      Call density (ns,ks,dens,p(1,2,i2),p(1,2,j1),'x')
      rk4 = rk4 + SUM_AmB(ns,ks,conv,dens,'x')

      Call mrk_qpqp(k)
      Call density (ns,ks,dens,p(1,2,i1),p(1,1,j2),'x')
      Call convol  (ns,ks,conv,dens,4,'s','s')
      Call density (ns,ks,dens,p(1,2,i2),p(1,1,j1),'x')
      rk4 = rk4 + SUM_AmB(ns,ks,conv,dens,'x')

      Call mrk_pqpq(k)
      Call density (ns,ks,dens,p(1,1,i1),p(1,2,j2),'x')
      Call convol  (ns,ks,conv,dens,4,'s','s')
      Call density (ns,ks,dens,p(1,1,i2),p(1,2,j1),'x')
      rk4 = rk4 + SUM_AmB(ns,ks,conv,dens,'x')

      End Function rk4


!======================================================================
      Real(8) Function hf_gk1 (i,j,k)
!======================================================================
!     Returns  Fk (i, j) base on the assembling two-electron
!     B-spline integrals (see module DBS_integral)
!----------------------------------------------------------------------
      Use DBS_grid, only: ns,ks         ! Rk(i,j,j,i)
!      Use df_orbitals
      Use DBS_orbitals_pq,  only: p => pq

      Implicit none
      Integer, intent(in) :: i,j,k
      Real(8), external :: SUM_AmB

      Real(8) :: xppi(ns,ns),xqqi(ns,ns),xpqi(ns,ns),xqpi(ns,ns), &
                 xppj(ns,ns),xqqj(ns,ns),xpqj(ns,ns),xqpj(ns,ns), &
                 xpppp(ns,ns),xqqqq(ns,ns),xqpqp(ns,ns),xpqpq(ns,ns)

      hf_gk1 = 0.d0

!     if(i.ne.igk.or.k.ne.kgk) then
        Call density (ns,ks,xppi,p(1,1,i),p(1,1,i),'x')
        Call mrk_pppp(k)
        Call convol (ns,ks,xpppp,xppi,4,'s','s')

        Call density (ns,ks,xqqi,p(1,2,i),p(1,2,i),'x')
        Call mrk_qqqq(k)
        Call convol (ns,ks,xqqqq,xqqi,4,'s','s')

        Call density (ns,ks,xpqi,p(1,1,i),p(1,2,i),'x')
        Call mrk_pqpq(k)
        Call convol (ns,ks,xpqpq,xpqi,4,'s','s')

        Call density (ns,ks,xqpi,p(1,2,i),p(1,1,i),'x')
        Call mrk_qpqp(k)
        Call convol (ns,ks,xqpqp,xqpi,4,'s','s')
!      end if

!      if(j.ne.jgk) then
       Call density (ns,ks,xppj,p(1,1,j),p(1,1,j),'x')
       Call density (ns,ks,xqqj,p(1,2,j),p(1,2,j),'x')
       Call density (ns,ks,xpqj,p(1,1,j),p(1,2,j),'x')
       Call density (ns,ks,xqpj,p(1,2,j),p(1,1,j),'x')
!      end if

      hf_gk1 =  SUM_AmB(ns,ks,xqqqq,xqqj,'x') + &
                SUM_AmB(ns,ks,xpppp,xppj,'x') + &
                SUM_AmB(ns,ks,xqpqp,xqpj,'x') + &
                SUM_AmB(ns,ks,xpqpq,xpqj,'x')



!      igk=i; jgk=j; kgk=k

      End Function hf_gk1




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
      Use DBS_integrals, only: rkb
      Use DBS_debug

      Implicit none
      Integer, intent(in) :: ns,ks,icase
      Character, intent(in) :: sym_i,sym_j
      Real(8), intent(in ) :: d(ns,*)
      Real(8), intent(out) :: a(ns,*)
      ! local variables
      Integer :: i,j, ip,jp, imin,imax, jmin,jmax, ii,jj
      Real(8) :: c,t1,t2, x(ns,ns)

      t1 = RRTC()

      if(icase.le.2) a(1:ns,1:ks)=0.d0
      if(icase.gt.2) a(1:ns,1:ns)=0.d0

      Select case(icase)
!----------------------------------------------------------------------
      Case(1)                                          !  I( . a ; . b)

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

        do ip=1,ks
        do i=1,ns-ip+1
           a(i,ip) = SUM(d(1:ns,1:ks)*rkb(i,1:ns,ip,1:ks))
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=0.d0
        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do ip=1,ks
        do i=ks+1-ip,ns;  c=0.d0
        do jp=1,ks
        do j=ks+1-jp,ns;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do ip=1,ks
        do i=ks+1-ip,ns;  c=0.d0
        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=0.d0
        do jp=1,ks
        do j=ks+1-jp,ns;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      end if
!----------------------------------------------------------------------
      Case(2)                                         !  I( a . ; b . )

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

        do ip=1,ks
        do i=1,ns-ip+1
           a(i,ip) = SUM(d(1:ns,1:ks)*rkb(1:ns,i,1:ks,ip))
        end do; end do

        do jp=1,ks
        do j=1,ns-jp+1;   c=0.d0
        do ip=1,ks
        do i=1,ns-ip+1;   c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=0.d0
        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns;  c=0.d0
        do ip=1,ks
        do i=ks+1-ip,ns;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=0.d0
        do ip=1,ks
        do i=ks+1-ip,ns;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns;  c=0.d0
        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      end if

!----------------------------------------------------------------------
      Case(3);  a(1:ns,1:ns) = 0.d0                    ! I( . a ; b . )

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

       do jp=1,ks
       do j =1,ns-jp+1; jj=j+jp-1
       do ip=1,ks
       do i =1,ns-ip+1; ii=i+ip-1
          c = rkb(i,j,ip,jp);      a( i,jj) = a( i,jj) + c*d(ii, j)
          if(ip.gt.1)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(jp.gt.1)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(ip.gt.1.and.jp.gt.1)  a(ii, j) = a(ii, j) + c*d( i,jj)
       end do;  end do
       end do;  end do



       Do ii = 1,ns
       Do jj = 1,ns
        x = 0.d0
        Do ip = 1,ks; i = ii-ip+ks; if(i.gt.ns) Cycle
        Do jp = 1,ks; j = jj-jp+ks; if(j.gt.ns) Cycle
         x(i,j) = rkb(i,j,ip,jp)
        End do; End do
        a(ii,jj) = SUM(d(1:ns,1:ns)*x(1:ns,1:ns))
       End do; End do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a( i,jj) = a( i,jj) + c*d(ii, j)
          if(ip.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(jp.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(ip.ne.ks.and.jp.ne.ks) a(ii, j) = a(ii, j) + c*d( i,jj)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

       do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax; jj=j+jp-ks
       do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
       do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);      a( i,jj) = a( i,jj) + c*d(ii, j)
       end do; end do
       end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks
        do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);      a( i,jj) = a( i,jj) + c*d(ii, j)
          if(jp.ne.ks)             a( i, j) = a( i, j) + c*d(ii,jj)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a( i,jj) = a( i,jj) + c*d(ii, j)
          if(ip.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
        end do; end do
        end do; end do

      end if
!----------------------------------------------------------------------
      Case(4);  a(1:ns,1:ns) = 0.d0                    ! I( a . ; . b )

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

       do jp=1,ks
       do j =1,ns-jp+1; jj=j+jp-1
       do ip=1,ks
       do i =1,ns-ip+1; ii=i+ip-1
          c = rkb(i,j,ip,jp);      a(ii, j) = a(ii, j) + c*d( i,jj)
          if(ip.gt.1)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(jp.gt.1)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(ip.gt.1.and.jp.gt.1)  a( i,jj) = a( i,jj) + c*d(ii, j)
       end do;  end do
       end do;  end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

       do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax; jj=j+jp-ks
       do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
       do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);      a(ii, j) = a(ii, j) + c*d( i,jj)
       end do; end do
       end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a(ii, j) = a(ii, j) + c*d( i,jj)
          if(ip.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(jp.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(ip.ne.ks.and.jp.ne.ks) a( i,jj) = a( i,jj) + c*d(ii, j)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks
        do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a(ii, j) = a(ii, j) + c*d( i,jj)
          if(jp.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a(ii, j) = a(ii, j) + c*d( i,jj)
          if(ip.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
        end do; end do
        end do; end do

      end if
!----------------------------------------------------------------------
      Case default

       Stop 'convol: unknown case'

      End Select

      t2 = RRTC()
      ic_convol = ic_convol + 1
      if(icase.le.2) time_convol=time_convol + (t2-t1)

      End Subroutine convol
