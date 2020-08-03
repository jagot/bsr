!=======================================================================
    Real(8) Function bvalu2 (t, bcoef, ns, ks, x, jderiv)
!=======================================================================
!   This routine is a modification of de Boor's bvalue routine, modified
!   to Return the value at right endpoint continuous from the left,
!   instead of 0. It assumes the usual knot multiplicity ks at the right
!   endpoint.
!   It calculates the value at  x of the jderiv-th derivative of spline
!   from b-representation. The spline is taken to be continuous from the
!   right.
!   SUBROUTINES contained: interv
!-----------------------------------------------------------------------
!   on entry
!   --------
!       t      knot sequence, of length  ns+ks, assumed nondecreasing.
!       bcoef  coefficient sequence, of length  ns .
!       ns     length of  bcoef, assumed positive.
!       ks     order of the spline .
!
!             . . . W A R N I N G . . .
!       The restriction  ks <= kmax (=15)  is imposed
!       arbitrarily by the parameter statement defining dimensions
!       for aj, dm, dp  below, but is  NEVER CHECKED.
!
!       x      the point at which to evaluate the spline.
!       jderiv integer giving the order of the derivative to be evaluated
!              ASSUMED to be zero or positive.
!
!   on exit
!   -------
!   bvalu2 - the value of the (jderiv)-th derivative of  f  at  x .
!
!   method
!   ------
!   the nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
!   cated with the aid of  INTERV. the  ks  b-coeffs of  f  relevant for
!   this interval are then obtained from  bcoef (or taken to be zero if
!   not explicitly available) and are then differenced  jderiv  times to
!   obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
!   precisely, with  j = jderiv, we have from x.(12) of the text that
!
!           (d**j)f  =  sum ( bcoef(.,j)*b(.,ks-j,t) )
!
!   where
!                  / bcoef(.),                     ,  j = 0
!                  /
!   bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
!                  / ----------------------------- ,  j > 0
!                  /    (t(.+ks-j) - t(.))/(ks-j)
!
!   then, we use repeatedly the fact that
!
!      sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
!   with
!                   (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
!      a(.,x)  =    ---------------------------------------
!                   (x - t(.))      + (t(.+m-1) - x)
!
!   to write  (d**j)f(x)  eventually as a linear combination of b-splines
!   of order  1 , and the coefficient for  b(i,1,t)(x)  must then
!   be the desired number  (d**j)f(x). (see x.(17)-(19) of text).
!-----------------------------------------------------------------------
    Implicit none
    Integer, Intent(in) :: ns, ks, jderiv
    Real(8), Intent(in) :: x, bcoef(ns), t(ns+ks)
    Real(8) :: aj(ks), dm(ks), dp(ks)
    Real(8) :: fkmj
    Integer :: i,mflag,km1,jcmin,imk,nmi,jcmax,jc,j,jj,kmj,ilo

    bvalu2 = 0.d0

    if (jderiv >= ks) Return

! ... find  i  such that 1 <= i < ns+ks  and  t(i) < t(i+1) and
! ... t(i) <= x < t(i+1) . if no such i can be found,  x  lies
! ... outside the support of  the spline  f  and  bvalu2 = 0.
! ... (the asymmetry in this choice of i makes f rightcontinuous)

    if( x /= t(ns+1) .OR. t(ns+1) /= t(ns+ks) ) then
      Call interv1(t,ns+ks,x, i, mflag)
      if (mflag /= 0) Return
    else
      i = ns
    end if

! ... if ks = 1 (and jderiv = 0), bvalu2 = bcoef(i).

    km1 = ks - 1
    if ( km1 <= 0 ) then
      bvalu2 = bcoef(i)
      Return
    end if

! ... store the ks b-spline coefficients relevant for the knot interval
! ... (t(i),t(i+1)) in aj(1),...,aj(ks); compute dm(j) = x - t(i+1-j),
! ... dp(j) = t(i+j) - x, j=1,...,ks-1. Set any of the aj not obtainable
! ... from input to zero. set any t.s not obtainable equal to t(1) or
! ... to t(ns+ks) appropriately.

    jcmin = 1
    imk = i - ks
    if (imk < 0) then
      jcmin = 1 - imk
      do j=1,i
        dm(j) = x - t(i+1-j)
      end do
      do j=i,km1
        aj(ks-j) = 0.
        dm(j) = dm(i)
      end do
    else
      do j=1,km1
        dm(j) = x - t(i+1-j)
      end do
    end if

    jcmax = ks
    nmi = ns - i
    if (nmi < 0) then
      jcmax = ks + nmi
      do j=1,jcmax
        dp(j) = t(i+j) - x
      end do
      do j=jcmax,km1
        aj(j+1) = 0.
        dp(j) = dp(jcmax)
      end do
    else
      do j=1,km1
        dp(j) = t(i+j) - x
      end do
    end if
      do jc=jcmin,jcmax
        aj(jc) = bcoef(imk + jc)
      end do

! ... difference the coefficients  jderiv  times.

    if (jderiv /= 0) then
      Do j=1,jderiv
        kmj = ks-j
        fkmj = kmj
        ilo = kmj
        Do jj=1,kmj
          aj(jj) = ((aj(jj+1) - aj(jj))/(dm(ilo) + dp(jj)))*fkmj
          ilo = ilo - 1
        End do
      End do
    end if

! ... compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
! ... given its relevant b-spline coeffs in aj(1),...,aj(ks-jderiv).

    if (jderiv /= km1) then
      Do j=jderiv+1,km1
        kmj = ks-j
        ilo = kmj
        Do jj=1,kmj
          aj(jj) = (aj(jj+1)*dm(ilo) + aj(jj)*dp(jj))/(dm(ilo)+dp(jj))
          ilo = ilo - 1
        End do
      End do
    end if
    bvalu2 = aj(1)

    END FUNCTION bvalu2

!=======================================================================
      Subroutine INTERV1 ( xt, lxt, x, left, mflag )
!=======================================================================
!
!     Computes  left = max( i ; 1 <= i <= lxt  .and.  xt(i) <= x )
!     which is the interval containing x .
!
!     A reformatted version of the de Boor routine
!
!      on entry
!      --------
!      xt  a Real sequence, of length lxt, assumed to be nondecreasing
!      lxt number of terms in the sequence xt .
!      x   the point whose location with respect to the sequence xt is
!          to be determined.
!
!      on exit
!      -------
!      left, mflag   integers, whose values are
!         1     -1   if               x .lt.  xt(1)
!         i      0   if   xt(i)  .le. x .lt. xt(i+1)
!        lxt     1   if  xt(lxt) .le. x
!
!      In particular,  mflag = 0 is the 'usual' case.  mflag .ne. 0
!      indicates that  x  lies outside the halfopen interval
!      xt(1) .le. y .lt. xt(lxt) . the asymmetric treatment of the
!      interval is due to the decision to make all pp functions cont-
!      inuous from the right. (left - ?)
!
!      ...  m e t h o d  ...
!
!  The program is designed to be efficient in the common situation where
!  it is called repeatedly, with  x  taken from an increasing or decrea-
!  sing sequence. this will happen, e.g., when a pp function is to be
!  graphed. the first guess for  left  is therefore taken to be the val-
!  ue Returned at the previous Call and stored in the  l o c a l  varia-
!  ble  ilo . a first check ascertains that  ilo .lt. lxt (this is nec-
!  essary since the present call may have nothing to do with the previ-
!  ous call). then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
!  ilo  and are done after just three comparisons.
!     otherwise, we repeatedly double the difference  istep = ihi - ilo
!  while also moving  ilo  and  ihi  in the direction of  x , until
!                      xt(ilo) .le. x .lt. xt(ihi) ,
!  after which we use bisection to get, in addition, ilo+1 = ihi .
!  left = ilo  is then Returned.
!-----------------------------------------------------------------------

      Integer :: left,lxt,mflag,ihi,ilo,istep,middle
      Real(8) :: x,xt(lxt)
      Data ilo /1/

      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt

   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100

!     ... now x .lt. xt(ilo) . decrease  ilo  to capture  x

      istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50

!     ... now x .ge. xt(ihi) . increase  ihi  to capture  x

   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt

!     ... now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval

   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100

!     note. it is assumed that middle = ilo in case ihi = ilo+1 .

      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50

!    ... set output and Return

   90 mflag = -1
      left = 1
                                        Return
  100 mflag = 0
      left = ilo
                                        Return
  110 mflag = 1
      left = lxt
                                        Return
      End Subroutine INTERV1

