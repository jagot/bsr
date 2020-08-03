!======================================================================
      Subroutine read_pqbs(nu)
!======================================================================
!     read B-spline w.f. from bsw-file (unit nu)
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,j,k,l,n,m,itype,nsw,ksw,mw,kp,kq
      Character(5) :: elw
      Integer, external :: Ifind_bsorb, Iadd_bsorb 
      Real(8) :: tt(ns+ks)

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

    1 read(nu,end=2) elw,mw
      Call EL_NLJK(elw,n,k,l,j,i)
      m = Ifind_bsorb(n,k,i,2) 
      mbs(m)=mw 
      pq(1:ns,1,m)=0.d0; read(nu) pq(1:mw,1,m)
      pq(1:ns,2,m)=0.d0; read(nu) pq(1:mw,2,m)
      bpq(:,1,m) = MATMUL(fpbs,pq(:,1,m))
      bpq(:,2,m) = MATMUL(fqbs,pq(:,2,m))
      go to 1
    2 Close(nu)

      End subroutine read_pqbs


!======================================================================
      Subroutine Read_dbsw(nu,mode,ishift)
!======================================================================
!     read B-spline r.w.f. from bsw-file (unit nu) 
!     mode = 0  - clear previuous w.f. in module DBS_orbitals_pq
!            1  - only read radial functions for existing orbitals
!            2  - read and add if needed
!     ishift -  shift the set indexes for new orbitals
!               (-1 -> nulify them)
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: nu, mode, ishift
      Integer :: i,j,k,l,n,m,itype,nsw,ksw,mw
      Integer, External :: Ifind_bsorb
      Character(5) :: elw
      Real(8) :: x

      if(mode.eq.0)  Call alloc_DBS_orbitals_pq(0,ns)

      rewind(nu)
      read(nu) itype,nsw,ksw
      if(ksw.ne.ks) Stop ' Read_dbsw:  ksw <> ks'
      if(nsw.ne.ns) Stop ' Read_dbsw:  nsw <> ns'
      if(grid_type.gt.0.and.itype.ne.grid_type) &
        Stop 'Stop in read_dbsw: another grid_type ?'

    1 read(nu,end=2) elw,mw
      Call EL_NLJK(elw,n,k,l,j,i); i=i+ishift; if(ishift.eq.-1) i=0
      m = Ifind_bsorb(n,k,i,0)
      if(mode.eq.1.and.m.eq.0) then; read(nu) x; read(nu) x; go to 1; end if
      m = Ifind_bsorb(n,k,i,2)
      pq(1:ns,1:2,m)=0.d0
      mbs(m)=mw
      read(nu) pq(1:mw,1,m)
      read(nu) pq(1:mw,2,m)

      bpq(1:ns,1,m) = MATMUL(fpbs, pq(:,1,m))
      bpq(1:ns,2,m) = MATMUL(fqbs, pq(:,2,m))

      go to 1
    2 Close(nu)

      End subroutine Read_dbsw
