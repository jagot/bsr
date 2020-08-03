!======================================================================
      Integer function Memory_splines()      Result(m)
!======================================================================
! ... provide required memory for splines arrays (in 4B words)
!----------------------------------------------------------------------
      Use spline_param, only: ks,ns,nv

      m = 0

! ... spline_densities:

      m = m + 2 * (ns*ks + 2*ns*ns) + 7

! ... spline_galerkin:

      m = m + 2 * (8*ns*ks + ns*ns)

! ... spline_grid:

      m = m + 2 * (ns+ks + 3*ks*ks*(nv+1) + 3*nv*ks)

! ... spline_hl:

      m = m + 2 * 2*ns*ks  + 1

! ... spline_integrals:

      m = m + 2 * ns*ns*(2*ks-1)*(2*ks-1) + 2

! ... spline_moments:

      m = m + 2 * (ks*ks*ks*ks*nv + 4*ks*ks*nv) + 2

! ... spline_integrals:

      m = m + 2 * (ns*ks*ns*ks + ns*ns*ks + ns*ks + ns*ns + ns) 

! ... spline_slater:

      m = m + 2 * (6*ns*ks + 2*ns + 2*ns*(3*ks-2)) + 6

! ... rk_triang:

      m = m + 2 * (5*ks + ks*ks + ks*ks*(ks+nv) + ks*ks*(ks+1)*(ks+1)/4 ) + 18

      End function Memory_splines
  



