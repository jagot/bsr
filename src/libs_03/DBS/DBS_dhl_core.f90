!=======================================================================
      Module dhl_core
!=======================================================================
! ... contain B-spline representation for one-electron Dirac hamiltonian
! ... including the interaction with closed core
!-----------------------------------------------------------------------
      Implicit none

      Real(8), allocatable :: hd(:,:,:)       !  < . | HD | . >

      Integer :: nkap = 0
      Integer, allocatable :: kap_list(:)

      Integer :: iv = 0
      Real(8), allocatable :: Lvec(:)

      Integer :: iv1 = 0, iv2 = 0
      Real(8) :: Lvalue

      Integer :: kcore_max = 7

      End Module dhl_core

!======================================================================
      Subroutine Alloc_dhl_core(ms,ncore)
!======================================================================
!     allocation arrays in dhl_core Module
!----------------------------------------------------------------------
      Use dhl_core
      Use DBS_orbitals_pq, only: nbf,kbs

      Implicit none
      Integer, intent(in) :: ms,ncore
      Integer :: i
      Integer, allocatable :: kp(:)
      Integer, external :: ipointer

      if(allocated(hd)) Deallocate(hd)
      if(allocated(kap_list)) Deallocate(kap_list);  nkap = 0
      if(allocated(Lvec)) Deallocate(Lvec); iv = 0

      if(ms.le.0) Return

      if(allocated(kp)) Deallocate(kp)
      Allocate(kp(nbf-ncore)); kp=0
      nkap=0
      Do i=ncore+1,nbf
       if(ipointer(nkap,kp,kbs(i)).ne.0) Cycle
       nkap=nkap+1; kp(nkap)=kbs(i)
      End do
      Allocate(kap_list(nkap))
      kap_list = kp(1:nkap)
      Deallocate(kp)

      Allocate(hd(ms,ms,nkap));  hd=0.d0
      Allocate(Lvec(ms)); iv=0

      End Subroutine Alloc_dhl_core


!======================================================================
      Subroutine Gen_dhl_core (ncore,mbreit,kappa)
!======================================================================
!     Generates HD-matrixes and integrals
!----------------------------------------------------------------------
      Use DBS_grid,        only: ms
      Use DBS_dhl_pq,      only: dhl
      Use dhl_core
      Use Timer

      Implicit none
      Integer, intent(in) :: ncore,mbreit,kappa
      Integer :: k,ik

      Real(8) :: t0, t1

      Call Alloc_dhl_core(ms,ncore)
      Call Alloc_bk4_data(0)

! ... set up hd matrix ...

      call CPU_time(t0)
      Do ik = 1,nkap;  k=kap_list(ik)
       if(kappa.ne.0.and.kappa.ne.k) Cycle
       call CPU_time(t1)
       Call MAT_dhl_pq(k)
       hd(1:ms,1:ms,ik) = dhl(1:ms,1:ms)
       call timed_section_now(t1, "MAT_dhl_pq")
      End do  ! over kappa
      call timed_section_now(t0, "Gen_dhl_core: setup hd matrix")

! ... core contribution:

      call CPU_time(t0)
      if(ncore.gt.0) then
         call CPU_time(t1)
         Call Add_int_core(ncore,mbreit,kappa)
         call timed_section_now(t0, "Gen_dhl_core: Add_int_core")
         call CPU_time(t1)
         Call Load_int_core
         call timed_section_now(t0, "Gen_dhl_core: Load_int_core")
      end if
      call timed_section_now(t0, "Gen_dhl_core: core contribution")

      End Subroutine Gen_dhl_core


!=====================================================================
      Subroutine Add_int_core(ncore,mbreit,kappa)
!=====================================================================
!     generate DHL-modification due to interaction with closed core
!---------------------------------------------------------------------
      Use dhl_core
      Use DBS_orbitals_pq, only: kbs
      Use bk4_data

      Implicit none
      Integer, intent(in) :: ncore,mbreit,kappa
      Integer :: i, kc,lc,jc, k,k1,k2, v, kk,ll,jj, ik, int
      Integer, external ::  l_kappa, j_kappa, itra
      Real(8) :: C,CC,CR, S(8)
      Real(8), external :: CJKJ, SMU

      Do ik=1,nkap; kk=kap_list(ik); ll=l_kappa(kk); jj=j_kappa(kk)
       if(kappa.ne.0.and.kappa.ne.kk) Cycle
      Do i=1,ncore; kc=kbs(i); lc=l_kappa(kc); jc=j_kappa(kc)

! ... direct interaction:

       C = dble(jc+1); k=0

       int=11; Call Add_bk4_data(k,int,i,ik,C)  !  Rk(Pi Pc; Pi Pc)
       int=21; Call Add_bk4_data(k,int,i,ik,C)  !  Rk(Qi Qc; Qi Qc)
       int=31; Call Add_bk4_data(k,int,i,ik,C)  !  Rk(Pi Qc; Pi Qc)
       int=41; Call Add_bk4_data(k,int,i,ik,C)  !  Rk(Qi Pc; Qi Pc)

! ... exchange interaction:

       k1=iabs(jc-jj)/2;  k2=iabs(jc+jj)/2
       Do k=k1,k2
        if(k.gt.kcore_max) Cycle

        C = -CJKJ(jc,k,jj)**2 / (jj+1)
        if(C.eq.0.d0) Cycle

        CR = C; if(mod(lc+k+ll,2).eq.1.or.itra(lc,k,ll).eq.0) CR = 0.d0

        if(CR.ne.0.d0) then
         int=13; Call Add_bk4_data(k,int,i,ik,C)   !  Rk(Pi Pc; Pc Pi)
         int=23; Call Add_bk4_data(k,int,i,ik,C)   !  Rk(Qi Qc; Qc Qi)
         int=33; Call Add_bk4_data(k,int,i,ik,C)   !  Rk(Pi Qc; Pc Qi)
         int=43; Call Add_bk4_data(k,int,i,ik,C)   !  Rk(Qi Pc; Qc Pi)
        end if

        if(mbreit.eq.0) Cycle

        Do v=k-1,k+1; if(v.lt.0) Cycle
         if(mod(lc+v+ll,2).ne.1) Cycle
         if(SMU(kk,kc,kc,kk,k,v,S).eq.0.d0) Cycle
         int=53; CC = C*(S(1)+S(4))
                 Call Add_bk4_data(v,int,i,ik,CC)  ! Sk(Pi Pc; Qc Qi)
         int=54; CC = C*(S(2)+S(3))
                 Call Add_bk4_data(v,int,i,ik,CC)  ! Sk(Pc Pi; Qi Qc)
         int=63; CC = C*(S(5)+S(6))
                 Call Add_bk4_data(v,int,i,ik,CC)  ! Sk(Pi Qc; Qc Pi)
         int=64; CC = C*(S(7)+S(8))
                 Call Add_bk4_data(v,int,i,ik,CC)  ! Sk(Pc Qi; Qj Pc)
        End do   ! over v
       End do    ! over k

      End do     ! over i
      End do     ! over ik

      End Subroutine Add_int_core


!=====================================================================
      Subroutine Load_int_core
!=====================================================================
!     evaluate the contribution to one-electron integrals the
!     interaction with core
!---------------------------------------------------------------------
      Use DBS_grid,         only: ns,ks,ms
      Use DBS_orbitals_pq,  only: pq
      Use DBS_convolutions
      Use dhl_core
      Use bk4_data

      Implicit none
      Character :: sym_i, sym_j, sym_d, sym_r
      Integer :: i, i1,i2, j1,j2, k,m, ik,io, iadd, int
      Real(8) :: C
      Real(8) ::  x(ns,ns),y(ns,ns),z(ns,ns)

      Do i=1,nbk;  k=kr1(i); io=kr3(i); ik=kr4(i)

       if(cbk(i).eq.0.d0) Cycle

       Select case(kr2(i))
        Case(11); Call mrk_pppp(k); sym_i='s'; sym_j='s'   !  Rk(Pi Pc; Pj Pc)
                  i1=1; i2=1; iadd=1; int=1
        Case(13); Call mrk_pppp(k); sym_i='s'; sym_j='s'   !  Rk(Pi Pc; Pc Pj)
                  i1=1; i2=1; iadd=1; int=3
        Case(21); Call mrk_qqqq(k); sym_i='s'; sym_j='s'   !  Rk(Qi Qc; Qj Qc)
                  i1=2; i2=2; iadd=2; int=1
        Case(23); Call mrk_qqqq(k); sym_i='s'; sym_j='s'   !  Rk(Qi Qc; Qc Qj)
                  i1=2; i2=2; iadd=2; int=3
        Case(31); Call mrk_pqpq(k); sym_i='s'; sym_j='s'   !  Rk(Pi Qc; Pj Qc)
                  i1=2; i2=2; iadd=1; int=1
        Case(33); Call mrk_pqpq(k); sym_i='s'; sym_j='s'   !  Rk(Pi Qc; Pc Qj)
                  i1=1; i2=2; iadd=3; int=3
        Case(41); Call mrk_qpqp(k); sym_i='s'; sym_j='s'   !  Rk(Qi Pc; Qj Pc)
                  i1=1; i2=1; iadd=2; int=1
        Case(43); Call mrk_qpqp(k); sym_i='s'; sym_j='s'   !  Rk(Qi Pc; Qc Pj)
                  i1=2; i2=1; iadd=4; int=3

        Case(53); Call msk_ppqq(k); sym_i='n'; sym_j='n'   !  Sk(Pi Pc; Qc Qj)
                  i1=2; i2=1; iadd=3; int=3
        Case(54); Call msk_ppqq(k); sym_i='n'; sym_j='n'   !  Sk(Pc Pi; Qj Qc)
                  i1=1; i2=2; iadd=4; int=4
        Case(63); Call msk_pqqp(k); sym_i='n'; sym_j='n'   !  Sk(Pi Qc; Qc Pj)
                  i1=2; i2=2; iadd=1; int=3
        Case(64); Call msk_pqqp(k); sym_i='n'; sym_j='n'   !  Sk(Pc Qi; Qj Pc)
                  i1=1; i2=1; iadd=2; int=4

        Case default;  Stop 'unknown integral in load_int_core'
       End Select

       Select case(int)
        Case(1); sym_d=sym_j; sym_r=sym_i       ! direct
        Case(2); sym_d=sym_i; sym_r=sym_j
        Case(3); sym_d='x';   sym_r='x'         ! exchange
        Case(4); sym_d='x';   sym_r='x'
       End Select

       m = 0
       if(i.eq.1) then
        m=1
       else
        if(kr1(i).ne.kr1(i-1)) m=1
        if(kr2(i).ne.kr2(i-1)) m=1
        if(kr3(i).ne.kr3(i-1)) m=1
        if(cbk(i-1).eq.0.d0)   m=1
       end if

       if(m.eq.1) then
        Call Density(ns,ks,x,pq(1,i1,io),pq(1,i2,io),sym_d)
        Call Convol(ns,ks,y,x,int,sym_i,sym_j)
        Call Full_mat_sym(ns,ks,y,z,sym_r)
       end if

       C = cbk(i)/2.d0;  x = C * z;  y = TRANSPOSE(x)

       Select case(iadd)
        Case(1); i1=1; i2=ns; j1=1; j2=ns
         hd(i1:i2,j1:j2,ik) = hd(i1:i2,j1:j2,ik) + x + y
        Case(2); i1=ns+1; i2=ms; j1=ns+1; j2=ms
         hd(i1:i2,j1:j2,ik) = hd(i1:i2,j1:j2,ik) + x + y
        Case(3); i1=1; i2=ns; j1=ns+1; j2=ms
         hd(i1:i2,j1:j2,ik) = hd(i1:i2,j1:j2,ik) + x
         hd(j1:j2,i1:i2,ik) = hd(j1:j2,i1:i2,ik) + y
        Case(4); i1=ns+1; i2=ms; j1=1; j2=ns
         hd(i1:i2,j1:j2,ik) = hd(i1:i2,j1:j2,ik) + x
         hd(j1:j2,i1:i2,ik) = hd(j1:j2,i1:i2,ik) + y
       End Select

      End do  ! over coeff.s

      End Subroutine Load_int_core


!====================================================================
      Subroutine Gen_Lvec(i,vv)
!====================================================================
!     Store hd(:,:) * P(i)  in Lvec
!====================================================================
      Use dhl_core
      Use DBS_grid, only: ns
      Use DBS_orbitals_pq, only: kbs

      Implicit none
      Integer, intent(in) :: i
      Real(8) :: vv(ns+ns)
      Integer :: k
      Integer, external :: Ipointer

      if(i.eq.iv) then; vv = Lvec; Return; end if

      k = Ipointer(nkap,kap_list,kbs(i))
      if(k.eq.0) Stop 'Gen_Lvec:  kappa = 0'
      Call Get_pv(i,vv,ns)
      Lvec = matmul(hd(:,:,k),vv(:))
      iv = i
      vv = Lvec

      End Subroutine Gen_Lvec


!====================================================================
      Real(8) Function Lval(i1,i2)
!====================================================================
!     Calculate and store P(i1) * hd(:,:) * P(i2)
!====================================================================
      Use dhl_core
      Use DBS_grid, only: ns
      Use DBS_orbitals_pq, only: kbs

      Implicit none
      Integer, intent(in) :: i1,i2
      Real(8) :: v1(ns+ns),v2(ns+ns)
      Integer :: k
      Integer, external :: Ipointer

      if(i1.eq.iv1.and.i2.eq.iv2) then; Lval = Lvalue; Return; end if

      if(kbs(i1).ne.kbs(i2)) then
       iv1=i1; iv2=i2; Lvalue=0.d0; Lval=0.d0; Return
      end if

      k = Ipointer(nkap,kap_list,kbs(i1))
      if(k.eq.0) Stop 'Lval:  kappa = 0'

      Call Get_pv(i1,v1,ns)
      v2 = matmul(hd(:,:,k),v1(:))
      Call Get_pv(i2,v1,ns)
      Lvalue =  SUM(v1(:)*v2(:))
      Lval = Lvalue; iv1=i1; iv2=i2

      End Function Lval

