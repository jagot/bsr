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

#ifdef DEBUG_SPEEDUPS
      Use Timer
#endif

      Implicit none
      Integer, intent(in)  :: nsw,ksw, nsv,ksv, ll,mm
      Real(8), intent(in)  :: tw(nsw+ksw), cw(nsw), &
                              bsp(nv+1,ks,ksv), fbs(ns,ns)
      Real(8), intent(out) :: cv(ns)
      ! local variables
      Integer :: i,j,ip,iv
      Real(8) :: a(ns,ns), ygw(nv,ks)
      Real(8), external :: bvalu2
      Integer :: ipiv(ns), info
#ifdef DEBUG_SPEEDUPS
      Real(8) :: cv_ref(ns)
      Real(8) :: discrepancy
      Real(8), parameter :: tolerance = sqrt(epsilon(1.d0))
#endif

! ... evaluate the function in gaussian points:

      call TimerStart('Convert_pq: bvalu2 omp')
      !$omp parallel do collapse(2) default(shared) private(i,j)
      Do i=1,nv
         Do j=1,ks
            ygw(i,j) = bvalu2(tw,cw,nsw,ksw,gr(i,j),0) * grw(i,j)
         End do
      End do
      !$omp end parallel do
      call TimerStop('Convert_pq: bvalu2 omp')

      call TimerStart('Convert_pq: bvalu2')
      Do i=1,nv
         Do j=1,ks
            ygw(i,j) = bvalu2(tw,cw,nsw,ksw,gr(i,j),0) * grw(i,j)
         End do
      End do
      call TimerStop('Convert_pq: bvalu2')

! ... form the vector of inner products of the radial function
! ... and the spline basis functions:

      call TimerStart('Convert_pq: inner products, old')
      cv_ref = 0.d0
      Do iv = 1,nv
         Do ip = 1,ksv
            i = iv+ip-1
            cv_ref(i) = cv_ref(i) + SUM(ygw(iv,:)*bsp(iv,:,ip))
         End do
      End do
      call TimerStop('Convert_pq: inner products, old')

      call TimerStart('Convert_pq: inner products, new')
      ! To parallelize the below loop, one would introduce a temporary
      ! matrix for cv and reduce over the extraneous dimension
      ! afterwards.
      cv = 0.d0
      Do iv = 1,nv
         cv(iv:iv+ksv-1) = cv(iv:iv+ksv-1) + matmul(transpose(bsp(iv,:,:)), ygw(iv,:))
      End do
      call TimerStop('Convert_pq: inner products, new')

      discrepancy = sum(abs(cv - cv_ref))
      write(*,'("Discrepancy: ",e26.16,", tolerance: ",e26.16)') discrepancy, tolerance
      if(discrepancy > tolerance) error stop "Discrepancy exceeds tolerance"

! ... B-spline overlap matrix:

      a(1:nsv-ll-mm,1:nsv-ll-mm) = fbs(ll+1:nsv-mm,ll+1:nsv-mm)

! ... solve the equation:  a cv = <B|cw>
#ifdef DEBUG_SPEEDUPS
      call TimerStart('Convert_pq: DGESV')
      cv_ref(:) = cv(:)
#endif
      Call dgesv(nsv-ll-mm, 1, a, ns, ipiv, cv(ll+1), ns, info)

#ifdef DEBUG_SPEEDUPS
      call TimerStop('Convert_pq: DGESV')
      write(*,*) "Running old implementation (Gauss-Jordan elimination) for comparison"

      a(1:nsv-ll-mm,1:nsv-ll-mm) = fbs(ll+1:nsv-mm,ll+1:nsv-mm)

      call TimerStart('Convert_pq: Gauss-Jordan elimination')
      Call gaussj (a,nsv-ll-mm,ns,cv_ref(ll+1),1,1)
      call TimerStop('Convert_pq: Gauss-Jordan elimination')

      discrepancy = sum(abs(cv - cv_ref))

      write(*,'("Discrepancy: ",e26.16,", tolerance: ",e26.16)') discrepancy, tolerance
      if(discrepancy > tolerance) error stop "Discrepancy exceeds tolerance"
#endif

      End Subroutine Convert_pq
