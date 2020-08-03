!=====================================================================
    SUBROUTINE mkgrid2(ns,ks,z,h,hmax,rmax,t)
!=====================================================================

!   sets up the knots for spline

      IMPLICIT NONE

      Integer(4), intent(in) :: ns,ks
      Real(8), intent(in) :: z,h,hmax,rmax
      Real(8), intent(out), dimension(ns+ks) :: t

! .. Local variables

      INTEGER:: n, i, m, ml, me, me1, me2, nv,nt
      REAL(KIND=8):: hp1, h2, tm, tx, hm

      ! .. determine ml, the number of distinct points from 0 to 1

       ml = 1.d0/h + 0.1
       hp1 = 1.d0 + h

      ! .. determine tmax

      tm = z*rmax
      hm = z*hmax

      ! .. determine final point of "exponential" grid
      ! .. me: number of points from 1 to (1+h)**me
      ! .. m:  number of points from (1+h)**me to tmax

      me1 = MAX(0.0d0, LOG(hm/h)/LOG(hp1)+1.d0)
      me2 = LOG(tm)/LOG(hp1)+1

      IF ( me2 <= me1 ) THEN
        me = me2
        m = 0
      ELSE
        me = me1
        tx = hp1**me
        h2 = h*tx/hp1
        m = NINT((tm-tx)/h2)
      END IF

      n = ml + me + m + ks -1
      if(n.ne.ns) Stop 'mkgrid2: ns <> n '
      nv = ns - (ks -1)
      nt = ns + ks

      ! .. establish the grid for z*r

      t = 0.d0

      DO i = ks+1, ks+ml
        t(i) = t(i-1) + h
      END DO

      DO i = ks+ml+1, ks+me+ml
        t(i) = t(i-1)*hp1
      END DO

      DO i = ks+me+ml+1, n+1
        t(i) = t(i-1) + h2
      END DO
      t(n+2:nt) = t(n+1)

      ! .. scale the values to the R variable

      t = t/z

    END SUBROUTINE mkgrid2

