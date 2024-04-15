!=========================================================================
      Subroutine bhwf(n,l,z,coef)
!=========================================================================
!     computes the spline expansion coefficients for hydrogenic radial
!     function for quantum numbers (n,l) and effective nuclear charge z.
!
!     SUBROUTINES called:
!         vhwf, vinty, and dpbsl (from LINPACK)
!
!     Calling sequence:
!
!              bhwf
!             ------
!            /  |   \
!         vhwf vinty dpbsl
!          ||        /   \
!         hnorm    daxpy ddot
!
!-------------------------------------------------------------------------
!
!     on entry
!     --------
!         n,l     orbital quantum numbers
!         z       effective nuclear charge
!
!     on exit
!     -------
!         coef    the spline coefficients defining the expansion of
!                 the radial function P(nl;r)
!
!-------------------------------------------------------------------------
      Use spline_param
      Use spline_grid
      Use spline_galerkin

      Implicit none
      Integer, intent(in) :: n, l
      Real(8), intent(IN) :: z
      Real(8), intent(out) :: coef(ns)

      ! .. Local variables

      Integer :: m, ierr
      Real(8) :: yr(nv,ks), bsl(ks,ns)

        ! .. obtain the values of the radial function
        !    at all the gaussian points
        ! .. and multiply by the gaussian weights

      Do m=1,ks

        CALL vhwf(n,l,z,nv,gr(1,m),yr(1,m))

        yr(:,m) = yr(:,m)*grw(:,m)

      End do

      ! .. form the vector of inner products of the radial function and the
      ! .. spline basis functions

      Call vinty(yr,coef)

      ! .. apply the boundary condition at the origin

      Call zfacsb(bsl,l)

      coef(1:l+1)=0.d0

      ! .. solve the system of equations bs coef = vy for coef

      Call dpbtrs ('U',ns,ks-1,1,bsl,ks,coef,ns,ierr)
      if(ierr.ne.0) Stop 'bhwf:  dpbtrs (LAPACK) failed'

      End Subroutine bhwf


!======================================================================
      SUBROUTINE zfacsb(bsl,l)
!======================================================================
!     Sets up the overlap bs which is a transpose of sb, <B_i,B_j>,
!     with the correct boundary condition at r=0 and factorizes bs.
!
!     SUBROUTINES called:  dpbfa (from LINPACK)
!----------------------------------------------------------------------
      Use spline_param
      Use spline_galerkin

      Implicit none
      Integer, intent(in) :: l
      Real(8), intent(out) :: bsl(ks,ns)

      ! .. local variables

      Integer :: i,j,ierr

      ! .. copy the array, converting to row oriented band storage mode

      bsl = TRANSPOSE(sb)

      ! .. apply boundary condition at r=0

      Do i = 1,l+1
       bsl(ks,i) = 1.d0
       Do j = 1,ks-1
         bsl(j,ks-j+i)=0.d0
       End do
      End do

!      CALL dpbfa(bsl,ks,ns,ks-1,ierr)
!      if (ierr /= 0 ) STOP 'facsb: dpbfa (LINPACK) failed'

      Call DPBTRF('U',ns,ks-1,bsl,ks,ierr)
      if (ierr.ne.0)  Stop 'facsb: dpbtrf (LAPACK) failed'


      End Subroutine zfacsb


!=========================================================================
      Subroutine vhwf(n,l,z,nr,r,vh)
!=========================================================================
!     This program returns the vector of values, vh(i), of a normalized
!     hydrogenic function with nuclear charge z and quantum numbers (n,l)
!     at values of the radius, r(i),  i=1,nr.
!
!     Calls:  hnorm
!-------------------------------------------------------------------------
!
!     on entry
!     --------
!         n     principal quantum number
!         l     angular quantum number
!         z     Nuclear charge
!         nr    number of points in the vector
!         r     values at which radial function is to be evaluated, in
!               increasing order
!
!     on exit
!     -------
!         vh    values of the hydrogenic radial functions
!
!-------------------------------------------------------------------------

      Implicit none
      Integer, intent(in) :: n, l, nr
      Real(8), intent(in) :: z
      Real(8), intent(in) :: r(nr)
      Real(8), intent(out) :: vh(nr)

      ! .. Local variables

      Integer :: k, i, nm
      Real(8) :: a, b, c, factor
      Real(8) :: w(nr), p(nr)
      Real(8), external :: hnorm1

      k = n-l-1
      if(k.lt.0) Stop ' vhwf:  n < l+1 '

      factor =  hnorm1(n,l,z)       ! .. gets the normalization 'factor'

      w = -2.D0*z*r/n  ! .. store the argument of the exponential factor

      ! .. to avoid underlow, determine point at which function will be zero
      ! .. Note that exp(-150) = 0.7175 10**-65

      nm = nr
      Do
        if (w(nm) < -150.D0 .AND. nm /= -1) then
          vh(nm) = 0.D0
          nm = nm-1
        else
          EXIT
        end if
      End do

      ! .. Initialize the recurrence relation for each value of r(i)

      p(1:nm) = 1.D0

      a = 1.D0
      b = k
      c = n+l

      ! .. Apply the recurrence relation when k is positive

      Do i = 1,k
        p(1:nm) = 1.D0+a/b*p(1:nm)/c*w(1:nm)
        a = a+1.D0
        b = b-1.D0
        c = c-1.D0
      End do

      ! .. Multiply the factor by the exponential

      vh(1:nm) = factor*p(1:nm)*EXP(w(1:nm)/2.D0)*(-w(1:nm))**(l+1)

      End Subroutine vhwf


!=====================================================================
      Real(8) Function hnorm1(n,l,Z)
!======================================================================
!     returns the value of the normalization constant for an
!     hydrogenic function with nuclear charge z, and orbital
!     quantum numbers (n,l)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,l
      Real(8), intent(in) :: z

      ! .. local variables

      Integer :: m, i
      Real(8) :: a, b, d, t

      m = l+l+1
      a = n+l
      b = m
      t = a
      d = b
      m = m-1

      if ( m > 0) then
        Do  i = 1,m
          a = a-1.d0
          b = b-1.d0
          t = t*a
          d = d*b
        End do
      end if

      hnorm1 = SQRT(z*t)/(n*d)

      End Function hnorm1

