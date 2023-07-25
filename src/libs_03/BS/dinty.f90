!======================================================================
      Subroutine DINTY (icase,ygr,ym)
!======================================================================
!
!    Computes the array elements   <B_i,  y(r) B _j>, if icase=0,
!                                  <B_i,  y(r) B'_j>, if icase=1,
!
!    in non-symmetric column storage mode
!
!    on entry
!    --------
!    ygr   array of values of a specific function  y(r) at the
!          gaussian points of each interval, weighted by the gaussian
!          weights
!
!    on exit
!    -------
!    ym    <B_i, y(r) B_j>  or  <B_i, y(r) B'_j>
!          in non-symmetric column storage mode
!--------------------------------------------------------------------
     Use spline_param
     Use spline_grid

     Implicit none
     Integer, intent(in) :: icase
     Real(8), intent(in) :: ygr(nv,ks)
     Real(8), intent(out) :: ym(ns,2*ks-1)

     ! .. local variables

     Integer :: ith, jth, i, irow, jcol

     ym = 0.d0                   ! clear ym array

     do i = 1,nv                 ! over intervals

      do ith = 1,ks             ! over left B-splines
        irow = i+ith-1
        do jth = 1,ks           ! over right B-splines
          jcol = jth-ith+ks

          if(icase.eq.0)  ym(irow,jcol) = ym(irow,jcol) + &
                SUM(ygr(i,:)*bsp(i,:,ith)*bsp (i,:,jth))

          if(icase.eq.1)  ym(irow,jcol) = ym(irow,jcol) + &
                SUM(ygr(i,:)*bsp(i,:,ith)*bspd(i,:,jth,1))

        end do
      end do

     end do

     End Subroutine DINTY
