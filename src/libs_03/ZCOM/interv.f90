!=======================================================================
      Subroutine INTERV ( xt, lxt, x, left, mflag )
!=======================================================================
!
!     Computes  left = max( i ; 1 <= i <= lxt  .and.  xt(i) <= x )
!     which is the interval containing x .
!
!     A reformatted version of the de Boor routine
!
!      on entry
!      --------
!
!      xt  a real sequence, of length lxt, assumed to be nondecreasing
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
!  ue returned at the previous call and stored in the  l o c a l  varia-
!  ble  ilo . a first check ascertains that  ilo .lt. lxt (this is nec-
!  essary since the present call may have nothing to do with the previ-
!  ous call). then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
!  ilo  and are done after just three comparisons.
!     otherwise, we repeatedly double the difference  istep = ihi - ilo
!  while also moving  ilo  and  ihi  in the direction of  x , until
!                      xt(ilo) .le. x .lt. xt(ihi) ,
!  after which we use bisection to get, in addition, ilo+1 = ihi .
!  left = ilo  is then returned.
!-----------------------------------------------------------------------

      integer left,lxt,mflag,   ihi,ilo,istep,middle
      double precision x,xt(lxt)
      data ilo /1/

      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt

   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100

!     ... now x .lt. xt(ilo) . decrease  ilo  to capture  x

   30 istep = 1
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

!    ... set output and return

   90 mflag = -1
      left = 1
                                        return
  100 mflag = 0
      left = ilo
                                        return
  110 mflag = 1
      left = lxt
                                        return
      End Subroutine INTERV
