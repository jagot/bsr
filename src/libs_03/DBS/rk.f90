!======================================================================
      Real(8) Function rk (f1,f2,f3,f4,k)
!======================================================================
!     Returns  rk - integral, base on the assembling two-electron
!     B-spline integrals (see module DBS_integral)
!----------------------------------------------------------------------
      Use DBS_grid,  only: ns,ks

      Implicit none
      Integer, intent(in) :: k
      Real(8), intent(in) :: f1(ns,2),f2(ns,2),f3(ns,2),f4(ns,2)
      Real(8) :: dens1(ns,ks),dens2(ns,ks),dens3(ns,ks),dens4(ns,ks), &
                 conv(ns,ks), S
      Real(8), external :: SUM_AmB

      rk = 0.d0

      Call density (ns,ks,dens1,f1(1,1),f3(1,1),'s')
      Call density (ns,ks,dens2,f2(1,1),f4(1,1),'s')
      Call density (ns,ks,dens3,f1(1,2),f3(1,2),'s')
      Call density (ns,ks,dens4,f2(1,2),f4(1,2),'s')

      Call mrk_pppp(k)
      Call convol  (ns,ks,conv,dens1,2,'s','s')
      S = SUM_AmB(ns,ks,conv,dens2,'s')
      rk = rk + S

      Call mrk_qqqq(k)
      Call convol  (ns,ks,conv,dens3,2,'s','s')
      S = SUM_AmB(ns,ks,conv,dens4,'s')
      rk = rk + S

      Call mrk_qpqp(k)
      Call convol  (ns,ks,conv,dens3,2,'s','s')
      S = SUM_AmB(ns,ks,conv,dens2,'s')
      rk = rk + S

      Call mrk_pqpq(k)
      Call convol  (ns,ks,conv,dens1,2,'s','s')
      S = SUM_AmB(ns,ks,conv,dens4,'s')
      rk = rk + S

      End Function rk

