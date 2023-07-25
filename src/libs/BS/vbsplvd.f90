
!=======================================================================
   SUBROUTINE vbsplvd(t, kg, ni, x, nderiv, dbiatx)
!=======================================================================
!
!  This routine calculates the values of the B-splines and their deriva-
!  tives, of order up to nderiv, that do not vanish at x(i), i=1,..ni
!  There are ks such B-splines at each point.
!
!  This routine is a vector version of bsplvd written by C. de Boor,
!  ``A Practical Guide to Splines".
!
!  subroutine contained: vbsplvb
!
!  calling sequence:
!       vbsplvd
!          ||
!       vbsplvb
!
!-----------------------------------------------------------------------
!  on entry
!  --------
!  t     the knot array, of length nt >=nv+2ks-1.  It is assumed
!        that t(i) < t(i+1) for each interval containing an x(i)
!        Division by zero will result otherwise (in vbsplvb).
!
!  kg    gives the beginning interval from which the B-splines are
!        evaluated at the Gaussian points.
!
!  ni    the number of intervals in which B-splines are to be evaluated
!        at all Gaussian points, its uplimit is nt.
!
!  x     the point array at which these values are sought,
!        one per interval, of length ni.
!
!  nderiv   an integer indicating that values of B-splines and their
!        derivatives up to but not including the  nderiv-th  are asked
!        for.
!
!  working area
!  ------------
!  w31   a three dimensional array, w31(i,j,m) (j=1,..,ks m=1,..,ks) con-
!        tains B-coeff.s of the derivatives of a certain order of the
!        ks B-splines of interest at point x(i)
!
!  w1,w2      one dimensional arrays
!
!  on return
!  ---------
!  dbiatx     a three dimensional array. its entry (i,j,m) contains
!        value of  (m-1)st  derivative of  (l-ks+j)-th B-spline of
!        order ks at point x(i) for knot sequence t, i=1..ni,
!        j=m..ks; m=1..nderiv;and l=kg..kg+ni-1
!
!  method
!  ------
!  values at x(i) of all the relevant B-splines of order ks,ks-1,...,
!  ks+1-nderiv  are generated via vbsplvb and stored temporarily
!  in dbiatx. then, the B-coeffs of the required derivatives of the
!  B-splines of interest are generated by differencing, each from the
!  preceding one of lower order, and combined with the values of B-
!  splines of corresponding order in dbiatx to produce the desired
!  values.
!----------------------------------------------------------------------

    USE spline_param

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kg, ni, nderiv
    REAL(KIND=8), DIMENSION(ns+ks), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(nv,ks,ks), INTENT(INOUT):: dbiatx
    REAL(KIND=8), DIMENSION(ni), INTENT(IN):: x

    ! local variables
    REAL(KIND=8):: fkpimm
    REAL(KIND=8), DIMENSION(ni):: w1,w2
    REAL(KIND=8), DIMENSION(ni,ks,ks):: w31
    REAL(KIND=8), DIMENSION(ni,ks) :: deltar, deltal
    INTEGER:: i, j, n, m, mhigh, kp1, jhigh, ideriv
    INTEGER:: ldummy, kpimm, jlow, jpimid, il

    mhigh = MAX(MIN(nderiv,ks),1)   !mhigh is usually equal to nderiv.
    kp1 = ks+1
    jhigh = kp1-mhigh
    CALL vbsplvb(kg,ni,x,jhigh,1,dbiatx)
    IF(mhigh == 1) RETURN

    ! ..the first row of dbiatx always contains the B-spline values
    ! ..for the current order. these are stored in row ks+1-current
    ! ..order before vbsplvb is called to put values for the next
    ! ..higher order on top of it. Vbsplvb only uses the first two dimensions

    ideriv = mhigh
    DO m = 2, mhigh
      jpimid = 1
      DO j = ideriv,ks
       dbiatx(1:ni,j,ideriv) = dbiatx(1:ni,jpimid,1)
       jpimid = jpimid+1
      END DO
      ideriv = ideriv-1
      jhigh = kp1-ideriv	
      CALL vbsplvb(kg,ni,x,jhigh,2,dbiatx)
    END DO

    ! at this point,  b(.,n-ks+i, ks+1-j)(x) is in dbiatx(.,i,j) for
    ! n=kg..kg+ni-1,i=j..ks,j=1..mhigh('='nderiv).in particular,the
    ! first row of  dbiatx  is already in final form. to obtain cor-
    ! responding derivatives of B-splines in subsequent rows, gene-
    ! rate their B-repr. by differencing, then evaluate at x(.).

    jlow = 1
    DO i = 1,ks
      w31(1:ni,jlow:ks,i) = 0.d0	
      jlow = i
      w31(1:ni,i,i) = 1.d0
    END DO

    ! at this point, w31(.,.,j) contains the B-coeffs for the j-th of the
    ! ks B-splines of interest here.

    DO m = 2,mhigh
      kpimm = kp1-m
      fkpimm = kpimm
      i = ks
      il = 0

      ! for j=1,...,ks, construct B-coeffs of  (m-1)st  derivative of
      ! B-splines from those for preceding derivative by differencing
      ! and store again in  w31(.,.,j). the fact that w31(i,j) = 0  for
      ! i < j is used.

      DO ldummy = 1, kpimm
       DO n = kg-il,ni+kg-il-1
        w1(n-kg+il+1) = fkpimm/(t(n+kpimm)-t(n))
       END DO

        ! the assumption that t(n) < t(n+1) makes denominator
        ! in w1(1..ni) nonzero.

        DO j = 1,i
         w31(1:ni,i,j) = (w31(1:ni,i,j)-w31(1:ni,i-1,j))*w1(1:ni)
        END DO
        il = il+1
        i = i-1
      END DO

      ! for i=1,...,ks, combine B-coeffs a(.,.,i) with B-spline values
      ! stored in dbiatx(.,.,m) to get value of (m-1)st  derivative of
      ! i-th B-spline (of interest here) at x(.), and store in
      ! dbiatx(.,i,m). storage of this value over the value of a B-spline
      ! of order m there is safe since the remaining B-spline derivat-
      ! ive of the same order do not use this value due to the fact
      ! that  a(.,j,i) = 0  for j .lt. i .

      DO i = 1,ks
        w2(1:ni) = 0.d0
        jlow = MAX(i,m)
        DO j = jlow,ks
          w2(1:ni) = w2(1:ni) + w31(1:ni,j,i)*dbiatx(1:ni,j,m)
        END DO
        dbiatx(1:ni,i,m) = w2(1:ni)
      END DO

    END DO

    CONTAINS


    !===================================================================
      SUBROUTINE vbsplvb(kg, ni, x, jhigh, index, biatx)
    !=====================================================================
    !  This routine calculates the values of all possibly nonzero B-splines
    !  at x(i) (i=1,..ni) of order
    !               jout=max(jhigh,(j+1)*(index-1))
    !  with knot sequence  t .
    !
    !  This routine is a vector version of bsplvb written by C. de Boor,
    !  "A Practical Guide to Splines", Chapter X, page 135
    !
    !  on entry
    !  --------
    !  t    -  knot sequence, of length nt=ns+ks, assumed to be nonde-
    !          creasing, that is t(i) <= t(i+1)
    !
    !  jhigh-  choose jhigh=1 to get the B-spline values directly
    !            by calling vbsplvb.
    !
    !  kg   -  gives the beginning interval from which the B-splines
    !           are to be evaluated at Gaussin points.
    !
    !  ni   -  the number of intervals in which B-splines are
    !            evaluated at all Gaussian points, its uplimit is nv.
    !
    !  x    -  the points at which the B-splines are to be evaluated,
    !            its length is ni,
    !
    !  index-  integers which determine the order  jout = max(jhigh,
    !            (j+1)*(index-1))  of the B-splines whose values at x(i)
    !            are to be returned.  index is used to avoid recalcula-
    !            tions when several columns of the triangular array of
    !            B-spline values are needed (e.g., in vbsplvd ).
    !            More precisely,
    !                     if index = 1 ,
    !            the calculation starts from scratch and the entire
    !            triangular array of B-spline values of orders
    !            1,2,...,jhigh is generated, order by order ,
    !            i.e., column by column .
    !                     if  index = 2 ,
    !            only the B-spline values of order  j+1, j+2, ..., jout
    !            are generated, the assumption being that  biatx,j,
    !            deltal,deltar are, on  entry, as they were on exit at the
    !            previous call. In particular, if  jhigh = 0, then
    !            jout = j+1, i.e., just the next column of B-spline
    !            values is generated.
    !
    !  working area
    !  ------------
    !  deltal, deltar: two dimensional arrays
    !  term, saved:    one dimensional arrays.
    !
    !  on return
    !  ---------
    !  biatx.....two dimensional array, with biatx(j-k+1,i)(j=k..ni)
    !        containing the value at x(j-k+1) of the polynomial of order
    !        jout which agrees with the B-spline b(j-jout+1,jout,t) on
    !        the interval (t(j),t(j+1)).
    !
    !  method
    !  ------
    !  The recurrence relation
    !
    !                       x - t(i)              t(i+j+1) - x
    !     b(i,j+1)(x)  =  -----------b(i,j)(x) + ---------------b(i+1,j)(x)
    !                     t(i+j)-t(i)            t(i+j+1)-t(i+1)
    !
    !  is used (repeatedly) to generate the (j+1)-vector  b(l-j,j+1)(x),
    !  ...,b(l,j+1)(x)  from the j-vector  b(l-j+1,j)(x),...,
    !  b(l,j)(x), storing the new values in  biatx  over the old. the
    !  facts that
    !            b(i,1) = 1  if  t(i) <= x < t(i+1)
    !  and that
    !            b(i,j)(x) = 0  unless  t(i) <= x < t(i+j)
    !  are used. the particular organization of the calculations follows al-
    !  gorithm  (8)  in chapter x of the text.
    !-----------------------------------------------------------------------

        USE spline_param
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(nv,ks), INTENT(INOUT):: biatx
        REAL(KIND=8), DIMENSION(nv), INTENT(IN):: x
        INTEGER, INTENT(IN):: kg, ni, jhigh, index

        ! .. Local variables
        INTEGER:: i, jp1, m
        INTEGER, SAVE:: j=1
        REAL(KIND=8), DIMENSION(ni) :: term, saved

        IF(index == 1) THEN
          j=1
          biatx(1:ni,1)=1.d0
          IF (j >= jhigh)  RETURN
        END IF

        DO
          jp1=j+1
          saved(1:ni)=0.d0

          DO i=1,ni
            deltar(i,j)=t(i+kg-1+j)-x(i)
            deltal(i,j)=x(i)-t(i+kg-j)
          END DO 	

          DO m=1,j
            DO i=1,ni
              term(i)=biatx(i,m)/(deltar(i,m)+deltal(i,jp1-m))
              biatx(i,m)=saved(i)+deltar(i,m)*term(i)
              saved(i)=deltal(i,jp1-m)*term(i)
            END DO
          END DO

          biatx(1:ni,jp1)=saved(1:ni)
          j=jp1
          IF (j >= jhigh) EXIT
        END DO
      END SUBROUTINE vbsplvb

  END SUBROUTINE vbsplvd