!=======================================================================
      Subroutine YVAL (id,iw,mm,yv,ygr)
!=======================================================================
!
!     This routine computes the values of  r^mm f(r), f'(r), f''(r)
!     at the gaussian points of each interval, where f(r) is defined
!     by the spline expansion vector yv.
!
!     on entry
!     --------
!     id      indivates derivatives of f(r)	(=0,1,2)
!     iw      if iw > 0, results are weighted by gaussian coef.s
!     mm      integer defining the power of r
!     yv      the spline expansion vector for the funciton f(r)
!
!     on exit
!     -------
!     ygr     array of values of r^mm f(r) at the gaussian points
!             of each interval
!-----------------------------------------------------------------------
      Use spline_param
      Use spline_grid

      Implicit none
      Integer, intent(in) :: id,iw,mm
      Real(8), intent(in) :: yv(ns)
      Real(8), intent(out) :: ygr(nv,ks)

      ! .. local variables

      Integer :: m, i, ith
      Real(8) :: gw(nv,ks)

      if(mm.eq. 0) then
         gw = 1.d0
      elseif(mm.eq. 1) then
         gw = gr
      elseif(mm.eq.-1) then
         gw = grm
      elseif(mm.gt. 1) then
         gw = gr**mm
      elseif(mm.lt.-1) then
         gw = grm**(-mm)
      end if

      if(iw.ne.0) gw = gw * grw

      ygr = 0.d0

      if(id.eq.0) then

      do m = 1,ks
        do i = 1,nv
          do ith = 1,ks
           ygr(i,m) = ygr(i,m) + yv(i+ith-1)*bsp(i,m,ith)*gw(i,m)
         end do
        end do
       end do

      elseif(id.eq.1.or.id.eq.2) then

       do m = 1,ks
        do i = 1,nv
         do ith = 1,ks
          ygr(i,m) = ygr(i,m) + yv(i+ith-1)*bspd(i,m,ith,id)*gw(i,m)
         end do
        end do
       end do

      end if

      End Subroutine YVAL
