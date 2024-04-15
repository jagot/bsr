!======================================================================
      Subroutine Convert_pq(nsw,ksw,tw,cw,nsv,ksv,cv,bsp,fbs,ll,mm)
!======================================================================
!     This programs converts B-spline orbital cw(1:nsw) given
!     on the grid tw to the B-spline orbital cv(1:nsv) with
!     grid defined in module DBS_grid
!     first ll and last mm splines are removed
!----------------------------------------------------------------------
      Use DBS_grid,   only: ns,ks,nv
      Use DBS_gauss,  only: gr, grw

      Implicit none
      Integer, intent(in)  :: nsw,ksw, nsv,ksv, ll,mm
      Real(8), intent(in)  :: tw(nsw+ksw), cw(nsw), &
                              bsp(nv+1,ks,ksv), fbs(ns,ns)
      Real(8), intent(out) :: cv(ns)
      ! local variables
      Integer :: i,j,ip,iv
      Real(8) :: a(ns,ns), ygw(nv,ks)
      Real(8), external :: bvalu2

! ... evaluate the function in gaussian points:

      Do i=1,nv; Do j=1,ks
       ygw(i,j) = bvalu2(tw,cw,nsw,ksw,gr(i,j),0) * grw(i,j)
      End do;  End do

! ... form the vector of inner products of the radial function
! ... and the spline basis functions:

      cv = 0.d0
      Do iv = 1,nv; Do ip = 1,ksv; i = iv+ip-1
       cv(i) = cv(i) + SUM(ygw(iv,:)*bsp(iv,:,ip))
      End do; End do

! ... B-spline overlap matrix:

      a(1:nsv-ll-mm,1:nsv-ll-mm) = fbs(ll+1:nsv-mm,ll+1:nsv-mm)

! ... solve the equation:  a cv = <B|cw>

      Call gaussj (a,nsv-ll-mm,ns,cv(ll+1),1,1)

      End Subroutine Convert_pq
