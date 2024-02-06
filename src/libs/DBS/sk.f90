!======================================================================
      Real(8) Function sk_ppqq (f1,f2,f3,f4,k)
!======================================================================
!     Returns  S^k(PPQQ) integral, based on the assembling two-electron
!     B-spline integrals (see module DBS_integrals).
!----------------------------------------------------------------------
      Use DBS_grid,          only: ns,ks

      Implicit none
      Integer, intent(in) :: k
      Real(8), intent(in) :: f1(ns,2),f2(ns,2),f3(ns,2),f4(ns,2)
      Real(8) :: di(ns,ns),dj(ns,ns),dc(ns,ns)
      Real(8), external :: SUM_AmB

      if(k.lt.0) Stop 'Stop in Sk_ppqq: k < 0'                 ! ???

      Call msk_ppqq(k)
      Call density (ns,ks,di,f1(1,1),f3(1,2),'n')
      Call convol  (ns,ks,dc,di,2,'n','n')
      Call density (ns,ks,dj,f2(1,1),f4(1,2),'n')
      sk_ppqq = SUM_AmB(ns,ks,dc,dj,'n')

      End Function sk_ppqq

!======================================================================
      Real(8) Function sk_pqqp (f1,f2,f3,f4,k)
!======================================================================
!     Returns  S^k (PQQP) integral, based on the assembling two-electron
!     B-spline integrals (see module DBS_integrals).
!----------------------------------------------------------------------
      Use DBS_grid,          only: ns,ks

      Implicit none
      Integer, intent(in) :: k
      Real(8), intent(in) :: f1(ns,2),f2(ns,2),f3(ns,2),f4(ns,2)
      Real(8) :: di(ns,ns),dj(ns,ns),dc(ns,ns)
      Real(8), external :: SUM_AmB

      if(k.lt.0) Stop 'Stop in sk_pqqp: k < 0'                   ! ???

      Call msk_pqqp(k)
      Call density (ns,ks,di,f1(1,1),f3(1,2),'n')
      Call convol  (ns,ks,dc,di,2,'n','n')
      Call density (ns,ks,dj,f2(1,2),f4(1,1),'n')
      sk_pqqp = SUM_AmB(ns,ks,dc,dj,'n')

      End Function sk_pqqp

