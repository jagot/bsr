!=====================================================================
      Subroutine ykf (j1,j2,k,yk)
!=====================================================================
!     Computes the spline solution of the differential equation
!
!           (d^2-k(k+1)/r^2)yk(r)= -(2k+1)(1/r) P_1(r)P_2(r)
!
!     on exit   (in module spline_slater)   ???
!     -------
!     yk    spline expansion of the solution
!     fyk   values at the gaussian points * (1/r)
!
!---------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use df_orbitals

      Implicit none
      Integer, intent(in) :: j1,j2,k
      Real(8) :: c, yp1(nv,ks),yq1(nv,ks),yp2(nv,ks),yq2(nv,ks), yr(nv,ks), &
                 yk(ns), a(ns,ns), a1(ns,ns), a2(ns,ns)
      REAL(8), external :: QUADR
      Integer :: i,j

! ... create the numerical values of the orbitals on the grid


      Call BVALUE_grid(ksp,p(1,1,j1),yp1,pbsp) 
      Call BVALUE_grid(ksq,p(1,2,j1),yq1,qbsp) 
      Call BVALUE_grid(ksp,p(1,1,j2),yp2,pbsp) 
      Call BVALUE_grid(ksq,p(1,2,j2),yq2,qbsp) 

! ... set up the array of yk(i) = INTEGRAL [(1/r)fc1(r)fc2(r)*B_i(r)]

      yr = grw*grm*(yp1*yp2+yq1*yq2)
 
      Call Vinty_pq(ks,yr,yk,bsp)

      c = - (2*k+1);   yk = c*yk

! ... set up the matrix (d^2-k(k+1)/r^2)

      Call ZINTYm (nv,ks,ks,ks,bsp,bspdd,grw,ns,a1)

      yr = grw*grm*grm

      Call ZINTYm (nv,ks,ks,ks,bsp,bsp,yr,ns,a2)

      c = -k*(k+1);  a = a1 + c*a2     

! ... boundary conditions:   ???

      C = 0.d0
      if(k.eq.0 ) C = QUADR(p(1,1,j1),p(1,1,j2),0)

! ... solve the matrix equation

      Do i=1,ns
       yk(i) = yk(i) - C * a(i,ns)
      End do 

      a1(1:ns-2,1:ns-2) = a(2:ns-1,2:ns-1)

      Call gaussj(a1,ns-2,ns,yk(2),1,1)

      yk(1)  = 0.d0
      yk(ns) = C

      End Subroutine ykf
