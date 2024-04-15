!======================================================================
      Real(8) Function rk_pq (i1,j1,i2,j2,k)
!======================================================================
!     Returns  rk_pq (i1, j1; i2, j2), base on the assembling two-electron
!     B-spline integrals (see module DBS_integral)
!----------------------------------------------------------------------
      Use DBS_grid,         only: ns,ks
      Use DBS_orbitals_pq,  only: p => pq

      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens1(ns,ks),dens2(ns,ks),dens3(ns,ks),dens4(ns,ks), &
                 conv(ns,ks)
      Real(8), external :: SUM_AmB

      rk_pq = 0.d0

      Call density (ns,ks,dens1,p(1,1,i1),p(1,1,i2),'s')
      Call density (ns,ks,dens2,p(1,1,j1),p(1,1,j2),'s')
      Call density (ns,ks,dens3,p(1,2,i1),p(1,2,i2),'s')
      Call density (ns,ks,dens4,p(1,2,j1),p(1,2,j2),'s')

      Call mrk_pppp(k)
      Call convol  (ns,ks,conv,dens1,2,'s','s')
      rk_pq = rk_pq + SUM_AmB(ns,ks,conv,dens2,'s')

      Call mrk_qqqq(k)
      Call convol  (ns,ks,conv,dens3,2,'s','s')
      rk_pq = rk_pq + SUM_AmB(ns,ks,conv,dens4,'s')

      Call mrk_qpqp(k)
      Call convol  (ns,ks,conv,dens3,2,'s','s')
      rk_pq = rk_pq + SUM_AmB(ns,ks,conv,dens2,'s')

      Call mrk_pqpq(k)
      Call convol  (ns,ks,conv,dens1,2,'s','s')
      rk_pq = rk_pq + SUM_AmB(ns,ks,conv,dens4,'s')

      End Function rk_pq


!======================================================================
      Real(8) Function rk_pppp (i1,j1,i2,j2,k)
!======================================================================
!     Returns  Rk (P_i1, P_j1; P_i2, P_j2)
!----------------------------------------------------------------------
      Use DBS_grid,         only: ns,ks
      Use DBS_orbitals_pq,  only: p => pq

      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens(ns,ks),conv(ns,ks)
      Real(8), external :: SUM_AmB

      Call mrk_pppp(k)
      Call density (ns,ks,dens,p(1,1,i1),p(1,1,i2),'s')
      Call convol  (ns,ks,conv,dens,2,'s','s')
      Call density (ns,ks,dens,p(1,1,j1),p(1,1,j2),'s')
      rk_pppp  = SUM_AmB(ns,ks,conv,dens,'s')

      End Function rk_pppp


!======================================================================
      Real(8) Function rk_qqqq (i1,j1,i2,j2,k)
!======================================================================
!     Returns  Rk (Q_i1, Q_j1; Q_i2, Q_j2)
!----------------------------------------------------------------------
      Use DBS_grid,         only: ns,ks
      Use DBS_orbitals_pq,  only: p => pq

      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens(ns,ks),conv(ns,ks)
      Real(8), external :: SUM_AmB

      Call mrk_qqqq(k)
      Call density (ns,ks,dens,p(1,2,i1),p(1,2,i2),'s')
      Call convol  (ns,ks,conv,dens,2,'s','s')
      Call density (ns,ks,dens,p(1,2,j1),p(1,2,j2),'s')
      rk_qqqq = SUM_AmB(ns,ks,conv,dens,'s')

      End Function rk_qqqq


!======================================================================
      Real(8) Function rk_qpqp (i1,j1,i2,j2,k)
!======================================================================
!     Returns  Rk (Q_i1, P_j1; Q_i2, P_j2)
!----------------------------------------------------------------------
      Use DBS_grid,         only: ns,ks
      Use DBS_orbitals_pq,  only: p => pq

      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens(ns,ks),conv(ns,ks)
      Real(8), external :: SUM_AmB

      Call mrk_qpqp(k)
      Call density (ns,ks,dens,p(1,2,i1),p(1,2,i2),'s')
      Call convol  (ns,ks,conv,dens,2,'s','s')
      Call density (ns,ks,dens,p(1,1,j1),p(1,1,j2),'s')
      rk_qpqp = SUM_AmB(ns,ks,conv,dens,'s')

      End Function rk_qpqp


!======================================================================
      Real(8) Function rk_pqpq (i1,j1,i2,j2,k)
!======================================================================
!     Returns  Rk (P_i1, Q_j1; P_i2, Q_j2)
!----------------------------------------------------------------------
      Use DBS_grid,      only: ns,ks
      Use DBS_orbitals_pq,  only: p  => pq

      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens(ns,ks),conv(ns,ks)
      Real(8), external :: SUM_AmB

      Call mrk_pqpq(k)
      Call density (ns,ks,dens,p(1,1,i1),p(1,1,i2),'s')
      Call convol  (ns,ks,conv,dens,2,'s','s')
      Call density (ns,ks,dens,p(1,2,j1),p(1,2,j2),'s')
      rk_pqpq = SUM_AmB(ns,ks,conv,dens,'s')

      End Function rk_pqpq

