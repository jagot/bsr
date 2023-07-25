!======================================================================
      Subroutine ZINTYN (ns,ks,nv,B1,B2,ygr,ym)
!======================================================================
!     Computes the array elements   <B1_i,  y(r) B2 _j>
!  
!     on entry
!     --------
!  
!     B1, B2 - two bases, spesified in gausian points
!  
!     ygr   array of values of a specific function  y(r) at the
!           gaussian points of each interval, weighted by the gaussian
!           weights
!  
!     on exit
!     -------
!     ym    <B1_i, y(r) B2_j>     in non-symmetric column storage mode
!--------------------------------------------------------------------
      Implicit none
   
      Integer, intent(in)  :: ns,ks,nv
      Real(8), intent(in)  :: B1(nv+1,ks,ks),B2(nv+1,ks,ks),ygr(nv,ks)
      Real(8), intent(out) :: ym(ns,ks+ks-1)
   
      ! .. local variables
   
      Integer :: ith, jth, i, irow, jcol
   
      ym(1:ns,1:ks+ks-1) = 0.d0       
      Do i = 1,nv                           ! over intervals
       Do ith = 1,ks; irow = i+ith-1        ! over left B1-splines
        Do jth = 1,ks; jcol = jth-ith+ks    ! over right B2-splines
         ym(irow,jcol) = ym(irow,jcol) + &
              SUM(ygr(i,:)*b1(i,:,ith)*b2(i,:,jth))
        End do
       End do
      End do
   
      End Subroutine ZINTYN
   