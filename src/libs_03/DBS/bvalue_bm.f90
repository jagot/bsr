!======================================================================
      Subroutine BVALUE_bm(ksm,a,y,BSP)
!======================================================================
!     Computes the function y(r) = SUM(i)  a_i * B_i(r)
!     in all gausian points + border values
!----------------------------------------------------------------------
      Use DBS_grid

      Implicit none
      Integer :: ksm,i,j,iv,ith,m
      Real(8), intent(in)  :: a(nv+ksm-1), bsp(nv+1,ks,ksm)
      Real(8), intent(out) :: y(nv*ks+2)

      y = 0.d0
      Do iv = 1,nv                 ! over intervals
       Do m = 1,ks                 ! over gausian points
        Do ith = 1,ksm             ! over B-splines in given interval
         i = iv+ith-1              ! B-spline index
         j = 1 + (iv-1)*ks + m     ! radial point index
         y(j) = y(j) + a(i)*bsp(iv,m,ith)
        End do
       End do
      End do

      m = nv*ks+2                  ! last point
      Do i=1,ksm
       y(m) = y(m) + bsp(nv+1,1,i) * a(ns-ks+i)
      End do

      m = 1                        ! first point
      Do i=1,ksm
       y(m) = y(m) + bsp(nv+1,2,i) * a(i)
      End do

      End Subroutine BVALUE_bm


!======================================================================
      Subroutine BVALUE_grid(ksm,a,y,BSP)
!======================================================================
!     Computes the function y(r) = SUM(i)  a_i * B_i(r)
!     in all gausian points 
!----------------------------------------------------------------------
      Use DBS_grid

      Implicit none
      Integer :: ksm,iv,ith,m
      Real(8), intent(in)  :: a(nv+ksm-1), bsp(nv+1,ks,ksm)
      Real(8), intent(out) :: y(nv,ks)

      y = 0.d0
      Do iv = 1,nv                 ! over intervals
       Do m = 1,ks                 ! over gausian points
        Do ith = 1,ksm             ! over B-splines in given interval
         y(iv,m) = y(iv,m) + a(iv+ith-1)*bsp(iv,m,ith)
        End do
       End do
      End do

      End Subroutine BVALUE_grid
