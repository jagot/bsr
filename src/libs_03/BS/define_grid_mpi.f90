!=========================================================================
      SUBROUTINE define_grid_mpi(z)
!=========================================================================
!     gets input data for the grid and sets up the spline knots 
!     Requires file  'knot.dat' with grid parameters 
! -------------------------------------------------------------------------
      Use MPI

      USE spline_param
      USE spline_grid

      IMPLICIT NONE

      Real(8) :: z
      Logical :: EX
      integer :: myid, ierr

! ... get input data for the grid

      Inquire(file=AF_grid,exist=EX)

      if(.not.EX) Call Stop_mpi(0,0, &
         'provide knot.dat file for splines parameters !')

      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Open(nug,file=AF_grid)
      Call Read_ipar_mpi(nug,'grid_type',grid_type)      
      Call Read_ipar_mpi(nug,'ns',ns)      
      Call Read_ipar_mpi(nug,'ks',ks)      
      Call Read_rpar_mpi(nug,'z',z)      
      Call Read_rpar_mpi(nug,'h',h)      
      Call Read_rpar_mpi(nug,'hmax',hmax)      
      Call Read_rpar_mpi(nug,'rmax',rmax)      

      Select case(grid_type)
       Case(0);  Call mkgrid (z) 
       Case(1);  Call mkgrid1(z) 
       Case default
       Call Stop_mpi(0,0,'Unknown grid_type ')
      End Select

      if(myid.eq.0) Call write_knot(z)

      Close(nug)

      END SUBROUTINE define_grid_mpi


!=========================================================================
      SUBROUTINE write_knot(z)
!=========================================================================

      USE spline_param
      USE spline_grid

      IMPLICIT NONE
      Double precision :: z

      rewind(nug)
      write(nug,'(a,i12,a)') 'ks = ',ks,' ==>  order  of splines (ks)'
      write(nug,'(a,i12,a)') 'ns = ',ns,' ==>  number of splines (ns)'
      write(nug,'(f12.5,a)') 'z = ',z,' ==>  nuclear charge (z)'
      write(nug,'(f12.5,a)') 'h = ',h,' ==>  step size from 0 to 1 (h for z*r, = (1/2)^n)'
      write(nug,'(f12.5,a)') 'hmax = ',hmax,' ==>  maximum step size (hmax for r)'
      write(nug,'(f12.5,a)') 'rmax = ',rmax,' ==>  maximum r (rmax)'
      write(nug,'(a)') '***'
      write(nug,'(5f12.5)') t(1:ns+ks)
      write(nug,'(a)') '***'
      write(nug,'(i12,a)') ml,' ==>  number of distinct knots from 0 to 1 (=1/h)'
      write(nug,'(i12,a)') me,' ==>  number of knots in the exponential region '
      write(nug,'(i12,a)') nv,' ==>  number of intervals '
      write(nug,'(f12.5,a)') (1/hmax)**2,' ==>  max k^2 (Ry) '


      END SUBROUTINE write_knot


!=====================================================================
      SUBROUTINE mkgrid(z)
!=====================================================================
!     sets up the knots for spline (original version)
!---------------------------------------------------------------------
      USE spline_param
      USE spline_grid

      IMPLICIT NONE

      INTEGER(4) :: n, i, m, me1, me2, nt
      REAL(8):: z,hp1, h2, tmax, tx, h0

      ! .. determine ml, the number of distinct points from 0 to 1

       h0 = h
       ml = 1.d0/h + 0.1
       h = 1.d0/ml
       hp1 = 1.d0 + h

      ! .. determine tmax

      tmax = z*rmax
      hmax = z*hmax

      ! .. determine final point of "exponential" grid
      ! .. me: number of points from 1 to (1+h)**me
      ! .. m:  number of points from (1+h)**me to tmax

      me1 = MAX(0.0d0, LOG(hmax/h)/LOG(hp1)+1.d0)
      me2 = LOG(tmax)/LOG(hp1)+1

      IF ( me2 <= me1 ) THEN
        me = me2
        m = 0
      ELSE
        me = me1
        tx = hp1**me
        h2 = h*tx/hp1
        m = NINT((tmax-tx)/h2)
      END IF

      n = ml + me + m + ks -1
      ns = n
      nv = ns - (ks -1)
      nt = ns + ks

      ! .. establish the grid for z*r

      ALLOCATE (t(nt)); t = 0.d0

      DO i = ks+1, ks+ml
        t(i) = t(i-1) + h0
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
      hmax = hmax/z

      END SUBROUTINE mkgrid


!=====================================================================
      SUBROUTINE mkgrid1(z)
!=====================================================================
!     sets up the knots for splines in semi-exponential grid:
!-----------------------------------------------------------------------
      USE spline_param
      USE spline_grid

      IMPLICIT NONE
      INTEGER(4):: i,nt
      Real(8) :: z 

      rmax = rmax * z
      hmax = hmax * z
      nt = rmax / h + 1 + ks
      if(Allocated(t)) Deallocate(t)
      ALLOCATE (t(nt))
      t = 0.d0
      ns = ks
	  
! ... determine ml, the number of intervals from 0 to 1
! ... first make sure that h = 1/n as recomended

      ml = 1.d0/h + 0.5;  h = 1.d0/ml
      Do i=1,ml 
       ns = ns + 1; t(ns) = t(ns-1) + h
      End do

! ... determine me, the number of intervals in "exponential" grid

      if(hmax.lt.h) hmax = h
      me = 0 
      Do 
       t(ns+1) = t(ns)*(1.d0+h)
       if(t(ns+1)-t(ns).gt.hmax) Exit
       ns = ns + 1
       if(t(ns).gt.rmax) Exit
       me = me + 1
      End do
        
! ... rest of interval with step = hmax

      if(t(ns).lt.rmax) then
       Do 
        t(ns+1) = t(ns) + hmax
        ns = ns + 1
        if(t(ns).ge.rmax) Exit
       End do
      end if

      if(t(ns).gt.rmax) then
       ns = ns - 1; me = me - 1
       t(ns+1:ns+ks) = rmax
       t(ns) = (t(ns+1)+t(ns-1))/2.d0
      end if

      nv = ns - ks + 1

! ... scale to the R variable

      t = t/z
      hmax = hmax / z
      rmax = rmax / z

      End Subroutine mkgrid1

