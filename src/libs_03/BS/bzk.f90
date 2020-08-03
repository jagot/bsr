!=======================================================================
      Subroutine BZK (j1,j2,k,kk)
!=======================================================================
!     Solves the following differential equation by the
!     spline-Galerkin method:
! 
!         (d + k/r) zk(r) = r^kk * P_j1(r) * P_j2(r),
! 
!     with boundary conditionas zk(0) = 0 if k >= 0,
!     and  dz/dr (0) = 0 if k = -1 (for Nk integrals)
! 
!     Calls:  FACDZK, YVAL, VINTY, dgbtrf (LAPACK) 
! 
!     on entry
!     --------
!     k      the order of zk
!     kk     defines right side of equation
!     j1,j2  orbital pointers
! 
!     on exit
!     -------
!     yk   B-spline reprezantation of Zk
!     fyk  Gauss-point reprezantation of Zk
!---------------------------------------------------------------------
      Use spline_param
      Use spline_grid
      Use spline_slater
      Use spline_orbitals, p => pbs
  
      Implicit none
      Integer, intent(in) :: j1,j2,k,kk
      Integer :: info
  
      if(kz.ne.k) Call FACDZK (ktx,k,dzk,ipvtz);  kz=k
  
!  ... create the numerical values of the orbitals on the grid
  
      if (iy1 .ne. j1) then
        Call YVAL (0,0,0,p(1,j1),fy1);  iy1=j1
      end if
  
      if (iy2 .ne. j2) then
        Call YVAL (0,0,0,p(1,j2),fy2);  iy2=j2
      end if
  
! ... set up the array of yk(i) = INTEGRAL 'r^kk*fc1(r)fc2(r)*B_i(r)'
 
      if(kk.eq.0) then
         fc = grw*fy1*fy2        ! Nk integrals
      elseif(kk.eq.1) then
         fc = gr*grw*fy1*fy2     ! Vk integral
      else
         Stop ' BZK: kk is wrong '
      end if
 
      Call Vinty(fc,yk)
 
! ... boundary conditions:
 
      yk(1)=0.d0;  if(k.eq.-1) yk(2)=0.d0   ! for Nk with k=-1
 
! ... solve the matrix equation
 
      Call DGBTRS ('N',ns,ks-1,ks-1,1,dzk,3*ks-2,ipvtz,yk,ns,info)
      if( info .ne. 0 ) Stop 'BZK: dgbtrs (LAPACK) failed'
 
! ... evaluates the function (1/r)zk(r) at all the gaussian points
 
      Call YVAL (0,0,-1,yk,fyk)
 
      End Subroutine BZK
 