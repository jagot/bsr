!=========================================================================
      Subroutine read_grid(z)
!=========================================================================
!     gets input data for the grid and sets up the spline knots 
!     Requires file  'knot.dat' with grid parameters 
! -------------------------------------------------------------------------
      Use spline_param
      Use spline_grid
      Implicit none
      Real(8) :: z

      Call Check_file(AF_grid)
      Open(nug,file=AF_grid)

      grid_type = 0
      Call Read_ipar(nug,'grid_type',grid_type)      

      Select case(grid_type)
       
       Case(0)

        rewind(nug)
        read(nug,*) ks
        read(nug,*) ns
        read(nug,*) z
        read(nug,*) h
        read(nug,*) hmax
        read(nug,*) rmax

        Call mkgrid(z)

       Case(1)

        Call Read_ipar(nug,'ks',ks)      
        Call Read_ipar(nug,'ns',ns)      
        Call Read_rpar(nug,'z' ,z)      
        Call Read_rpar(nug,'h' ,h)      
        Call Read_rpar(nug,'hmax',hmax)      
        Call Read_rpar(nug,'rmax',rmax)      

        Call mkgrid1(z)

      End Select

      Close(nug)

      End Subroutine read_grid



!=========================================================================
      Subroutine define_grid(z)
!=========================================================================
!     gets input data for the grid and sets up the spline knots 
!     Requires file  'knot.dat' with grid parameters 
! -------------------------------------------------------------------------
      Use spline_param
      Use spline_grid
      Implicit none
      Real(8) :: z
      Logical :: EX

! ... get input data for the grid

      Inquire(file=AF_grid,exist=EX)

      if(.not.EX) then
       Call mkgrid1(z)
       Open(nug,file=AF_grid)
       Call Print_grid_01(z)
       Stop 'check knot.dat file for splines parameters ...'
      End if 

      Open(nug,file=AF_grid); Rewind(nug)

      grid_type = 0
      Call Read_ipar(nug,'grid_type',grid_type)      

      Select case(grid_type)
       
       Case(0);  Call define_grid_00(z) 
       Case(1);  Call define_grid_01(z) 
       Case default; Stop ' Unknown grid_type '

      End Select

      Close(nug)

      End Subroutine define_grid


!=========================================================================
      Subroutine define_grid_00(z)
!=========================================================================
      Use spline_param
      Use spline_grid

      Implicit none

      Real(8) :: z

      rewind(nug)
      read(nug,*) ks
      read(nug,*) ns
      read(nug,*) z
      read(nug,*) h
      read(nug,*) hmax
      read(nug,*) rmax

! ... set up the knots for spline

      Call mkgrid(z)

! ... print grid in the file 'knot.dat'

      rewind(nug)
      write(nug,'(i12,a)') ks,' ==>  order  of splines (ks)'
      write(nug,'(i12,a)') ns,' ==>  number of splines (ns)'
      write(nug,'(f12.5,a)') z,' ==>  nuclear charge (z)'
      write(nug,'(f12.5,a)') h,' ==>  step size from 0 to 1 (h for z*r, = (1/2)^n)'
      write(nug,'(f12.5,a)') hmax,' ==>  maximum step size (hmax for r)'
      write(nug,'(f12.5,a)') rmax,' ==>  maximum r (rmax)'
      write(nug,'(a)') '***'
      write(nug,'(5f12.5)') t(1:ns+ks)
      write(nug,'(a)') '***'
      write(nug,'(i12,a)') ml,' ==>  number of distinct knots from 0 to 1 (=1/h)'
      write(nug,'(i12,a)') me,' ==>  number of knots in the exponential region '
      write(nug,'(i12,a)') nv,' ==>  number of intervals '
      write(nug,'(f12.5,a)') (1/hmax)**2*3,' ==>  max k^2 (Ry) '


      End Subroutine define_grid_00


!=====================================================================
      Subroutine mkgrid(z)
!=====================================================================
!     sets up the knots for spline (original version)
!---------------------------------------------------------------------
      Use spline_param
      Use spline_grid

      Implicit none
      Integer :: n, i, m, me1, me2, nt
      Real(8) :: z,hp1, h2, tmax, tx, h0

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
      End IF

      n = ml + me + m + ks -1
      ns = n
      nv = ns - (ks -1)
      nt = ns + ks

      ! .. establish the grid for z*r

      Allocate (t(nt)); t = 0.d0

      Do i = ks+1, ks+ml
        t(i) = t(i-1) + h0
      End do

      Do i = ks+ml+1, ks+me+ml
        t(i) = t(i-1)*hp1
      End DO

      Do i = ks+me+ml+1, n+1
        t(i) = t(i-1) + h2
      End do
      t(n+2:nt) = t(n+1)

      ! .. scale the values to the R variable

      t = t/z
      hmax = hmax/z

      End Subroutine mkgrid


!=========================================================================
      Subroutine define_grid_01(z)
!=========================================================================
      Use spline_param
      Use spline_grid

      Implicit none
      Real(8) :: z

      Call Read_ipar(nug,'ks',ks)      
      Call Read_ipar(nug,'ns',ns)      
      Call Read_rpar(nug,'z' ,z)      
      Call Read_rpar(nug,'h' ,h)      
      Call Read_rpar(nug,'hmax',hmax)      
      Call Read_rpar(nug,'rmax',rmax)      

! ... set up the knots for spline

      CALL mkgrid1(z)

! ... print grid in the file 'knot.dat'

      Call Print_grid_01(z)

      End Subroutine define_grid_01


!=========================================================================
      Subroutine Print_grid_01(z)
!=========================================================================
      Use spline_param
      Use spline_grid

      Implicit none
      Real(8) :: z

      rewind(nug)
      write(nug,'(a,i3,T20,a)') 'grid_type  = ',grid_type,' !  type of grid '
      write(nug,*)
      write(nug,'(a,i10,T20,a)')   'ks   =',ks,' ! order of splines '
      write(nug,'(a,i10,T20,a)')   'ns   =',ns,' ! number of splines '
      write(nug,'(a,f10.5,T20,a)') 'z    =',z, ' ! nuclear charge'
      write(nug,'(a,f10.5,T20,a)') 'h    =',h, ' ! step size from 0 to 1, for z*r'
      write(nug,'(a,f10.5,T20,a)') 'hmax =',hmax,' ! maximum step size for r '
      write(nug,'(a,f10.5,T20,a)') 'rmax =',rmax,' ! maximum r '
      write(nug,'(a)') '***'
      write(nug,'(5f12.5)') t(1:ns+ks)
      write(nug,'(a)') '***'
      write(nug,'(i12,a)') ml,' ==>  number of distinct knots from 0 to 1 (=1/h)'
      write(nug,'(i12,a)') me,' ==>  number of knots in the exponential region '
      write(nug,'(i12,a)') nv,' ==>  number of intervals '
      write(nug,'(f12.5,a)') (1/hmax)**2*3,' ==>  max k^2 (Ry) '

      End Subroutine Print_grid_01


!=====================================================================
      Subroutine mkgrid1(z)
!=====================================================================
!     sets up the knots for splines in semi-exponential grid:
!-----------------------------------------------------------------------
      Use spline_param
      Use spline_grid

      Implicit none
      Integer:: i,nt
      Real(8) :: z 

      rmax = rmax * z
      hmax = hmax * z
      nt = (rmax / h + 1 + ks)*2
      if(Allocated(t)) Deallocate(t)
      Allocate (t(nt))
      t = 0.d0
      ns = ks
	  
! ... determine ml, the number of intervals from 0 to 1
! ... first make sure that h = 1/n as recomEnded

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
      End if

      if(t(ns).gt.rmax) then
       ns = ns - 1; me = me - 1
       t(ns+1:ns+ks) = rmax
       t(ns) = (t(ns+1)+t(ns-1))/2.d0
      End if

      nv = ns - ks + 1

! ... scale to the R variable

      t = t/z
      hmax = hmax / z
      rmax = rmax / z

      End Subroutine mkgrid1




