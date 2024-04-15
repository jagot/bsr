!======================================================================
      Module DBS_grid
!======================================================================
!     the B-spline knot parameters
!----------------------------------------------------------------------
      Implicit none

      Integer :: grid_type = 1    ! type of grid

      Integer :: ks = 9           ! order of B-splines
      Integer :: ns = 0           ! number of splines
      Integer :: ms = 0           ! number of splines in (P,Q) basis
      Integer :: nv = 0           ! number of intervals ( = ns-ks+1 )

      Integer :: ml = 1           ! number of initial equal intervals
      Integer :: me = 0           ! number of intervals in the exponential region

      Real(8) :: hi   =  0.25d0   ! initial step parameter
      Real(8) :: he   =  0.25d0   ! exponential step factor
      Real(8) :: hmax =  1.d0     ! maximum step, t(ns+1) - t(ns)
      Real(8) :: rmax = 50.d0     ! input border radius
      Real(8) :: tmax =  0.d0     ! real border radius, t(ns+1)

! ... additional bases:

      Integer :: ksp=8, nsp=0     ! B-splines for p-functions
      Integer :: ksq=9, nsq=0     ! B-splines for q-functions

! ... knot sequence, t(1:ns+ks)

      Real(8), allocatable :: t(:)

! ... where to look for the grid:

      Character(80) :: AF_grid = 'knot.dat'

      End Module DBS_grid

!======================================================================
      SUBROUTINE def_grid (knot,name,z,atw)
!======================================================================
!     get input data for the grid and sets up the spline knots
!
!     knot -  if not empty, it defines file with grid parameters
!     name -  name of case; if is not empty, will be created name.knot,
!             file with grid parameters for the given case
!     z    -  nuclear charge
!     atw  -  atomic weight
!
!     if "knot" file is empty, program check file "name.knot";
!     if both files are absent, program creats grid file for given z;
!     if z=0, just example file created for z=1 and execution stops
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_nuclear

      Implicit none
      Character(*) :: knot, name
      Real(8) :: z,atw
      Integer :: nug=0
      Logical :: EX

! ... define the file name for input grid parameters:

      if(len_trim(knot).ne.0) then
       AF_grid = knot
       Inquire(file=AF_grid,exist=EX)
       if(EX) Call Find_free_unit(nug)
       if(INDEX(knot,'.bsw').eq.0) then
        if(nug.gt.0) Open(nug,file=AF_grid)
       else
        if(nug.gt.0) Open(nug,file=AF_grid,form='UNFORMATTED')
        grid_type = -2
       end if
      end if

      if(nug.eq.0.and.len_trim(name).ne.0) then
       AF_grid = trim(name)//'.knot'
       Call Find_free_unit(nug)
       Open(nug,file=AF_grid)
      end if

! ... create knot.dat as example:

      if(z.eq.0.and.nug.eq.0) then
       z = 1.d0; atw = 1.d0
       Call Read_nuclear(0,z,atw)
       Call mkgrid_01
       Call Write_knotdat
       Stop 'check knot.dat file first'
      end if

! ... check nuclear parameters first:

      if(grid_type.ne.-2) then
       Call Read_nuclear(nug,z,atw)
      else
       Call Read_nuclear(0,z,atw)
      end if

! ... read grid parameters from the file if any:

      if(nug.gt.0) then
       Call Read_ipar(nug,'grid_type',grid_type)
       Call Read_ipar(nug,'ns',ns)
       Call Read_ipar(nug,'ks',ks)
       Call Read_ipar(nug,'nv',nv)
       Call Read_ipar(nug,'ml',ml)
       Call Read_rpar(nug,'hi',hi)
       Call Read_rpar(nug,'he',he)
       Call Read_rpar(nug,'hmax',hmax)
       Call Read_rpar(nug,'rmax',rmax)
       Call Read_ipar(nug,'ksp',ksp)
       Call Read_ipar(nug,'ksq',ksq)
      end if

! ... read grid parameters from command line if any:

      Call Read_iarg('grid_type',grid_type)
      Call Read_iarg('ns',ns)
      Call Read_iarg('ks',ks)
      Call Read_iarg('nv',nv)
      Call Read_iarg('ml',ml)
      Call Read_rarg('hi',hi)
      Call Read_rarg('he',he)
      Call Read_rarg('hmax',hmax)
      Call Read_rarg('rmax',rmax)
      Call Read_iarg('ksp',ksp)
      Call Read_iarg('ksq',ksq)

      if(ksp.le.0.or.ksp.gt.ks) ksp=ks-1
      if(ksq.le.0.or.ksq.gt.ks) ksq=ks

! ... create the knot grid:

      Select case(grid_type)
       Case(1);      Call mkgrid_01
       Case(2);      Call mkgrid_02
       Case(-1);     Call read_grid(nug)
       Case(-2);     Call read_grid_bsw(nug)
       Case default; Stop 'Unknown grid_type'
      End Select

      Close(nug)

      ms = ns + ns
      nv = ns - ks + 1
      nsp=nv+ksp-1
      nsq=nv+ksq-1
      tmax = t(ns+1)

! ... record the grid parameters:

      if(grid_type.gt.0) Call Write_knotdat

      End Subroutine def_grid


!======================================================================
      Subroutine read_knot_dat
!======================================================================
      Character(40) :: knot = 'knot.dat', name = ' '
      Real(8) :: z = 0.d0, awt = 0.d0

      Call def_grid(knot,name,z,awt)

      End  Subroutine read_knot_dat


!======================================================================
      Subroutine read_grid (nu)
!======================================================================
!     read t-sequence from knot.dat;
!     it is option to work with old knot-sequences
!----------------------------------------------------------------------
      Use DBS_grid

      Implicit none
      Integer :: nu, i,j
      Integer, external :: Ifind_position

      if(ns.le.0)  Stop 'Stop in read_grid:  ns = 0'
      if(allocated(t)) Deallocate(t); Allocate(t(ns+ks))
      if(Ifind_position(nu,'grid points:').eq.0) &
       Stop 'Stop in read_grid:  can not find grid points'

      read(nu,*)
      read(nu,*)
      Do i=1,ns+ks; read(nu,*) j,t(j); End do

      End Subroutine read_grid


!======================================================================
      Subroutine read_grid_bsw (nu)
!======================================================================
!     read t-sequence from name.bsw unformated files;
!     it is option to work with old knot-sequences
!----------------------------------------------------------------------
      Use DBS_grid

      Implicit none
      Integer :: nu, i

      rewind(nu)
      read(nu) i,ns,ks
      if(allocated(t)) Deallocate(t); Allocate(t(ns+ks))
      rewind(nu)
      read(nu) i,ns,ks,t,ksp,ksq

      End Subroutine read_grid_bsw



!======================================================================
      Subroutine Write_knotdat
!======================================================================
!     create example of the knot.dat file
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_nuclear

      Implicit none
      Integer :: i, nug

      Call Find_free_unit(nug)
      Open(nug,file=AF_grid)
      rewind(nug)

      write(nug,'(a,i3,6x)') 'grid_type = ',grid_type
      write(nug,*)
      write(nug,'(a,i5,15x,a)') 'ks =',ks,'!  order  of splines'
      write(nug,'(a,i5,15x,a)') 'ns =',ns,'!  number of splines'
      write(nug,'(a,i5,15x,a)') 'nv =',nv,'!  number of intervals'
      write(nug,'(a,i5,15x,a)') 'ml =',ml,'!  number of initial equal-spaced intervals'
      write(nug,*)
      write(nug,'(a,i5,15x,a)') 'ksp=',ksp
      write(nug,'(a,i5,15x,a)') 'ksq=',ksq
      write(nug,*)
      write(nug,'(a,f16.8,2x,a)') 'hi =   ',hi,  '!  initial step for r*zg '
      write(nug,'(a,f16.8,2x,a)') 'he =   ',he,  '!  exponetial step size '
      write(nug,'(a,f16.8,2x,a)') 'hmax = ',hmax,'!  maximum step size for r, given'
      write(nug,'(a,f16.8,2x,a)') 'rmax = ',rmax,'!  maximum radius, given'
      write(nug,'(a,f16.8,2x,a)') 'tmax = ',tmax,'!  maximum radius, actual'
      write(nug,*)
      write(nug,'(a,f8.4)') 'atomic_number = ',atomic_number
      write(nug,'(a,f8.4)') 'atomic_weight = ',atomic_weight
      write(nug,'(a,f8.4)') 'rrms          = ',rrms
      write(nug,*)
      write(nug,'(a,a)') 'nuclear =  ',nuclear

      if(nuclear.eq.'Fermi') then
       write(nug,*)
       write(nug,'(a,d24.16,a)') 'c_fermi =   ',c_fermi,' !  in a.u.'
       write(nug,'(a,d24.16,a)') 'a_fermi =   ',a_fermi,' !  in a.u.'
      end if

      if(nuclear.eq.'uniform') then
       write(nug,*)
       write(nug,'(a,E16.8,a)') 'r_uniform = ',r_uniform
      end if

      write(nug,*)
      write(nug,'(a)') 'grid points:'
      write(nug,*)
      if(allocated(t)) then
       write(nug,'(i5,d24.16)') (i,t(i),i=1,ns+ks)
      End if
      write(nug,'(a)') '***'
      write(nug,*)
      write(nug,'(a)') 'possible grid-type:'
      write(nug,*)
      write(nug,'(a)') ' 1  -  semi-exponential, ns is derived from the parameters  '
      write(nug,'(a)') ' 2  -  exponential, for given ns and rmax '
      write(nug,'(a)') '-1  -  read from the given file '
      write(nug,'(a)') '-2  -  read from the bsw-file '

      Close(nug)

      End Subroutine Write_knotdat


!======================================================================
      Subroutine mkgrid_01
!======================================================================
!     sets up the knots for splines in semi-exponential grid:
!      1. first ml equal intervals "hi"  (this part depends on nuclear)
!      2. then exponantial grid with interval increasing as (1+he)
!      3. then again equal intervals "hmax" if any
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_nuclear, z => atomic_number

      Implicit none
      Integer :: i,ne
      Real(8) :: ti,tj,ts,ht,hs

! ... initial point:

      if(nuclear.eq.'point') then
       if(z.le.0.d0) Stop 'mkgrid_01: Z <= 0'
       hs=hi/z
      elseif(nuclear.eq.'Fermi') then
       ts = 4.d0 * log(3.d0) * a_fermi
       hs = ts*hi
      else
       hs = r_uniform*hi
      End if

      ti = hs*ml; ns = ks + ml

! ... "exponential" part of grid:

      ht = 1.d0 + he
      Do
       tj = ti*ht
       if(tj-ti.gt.hmax) Exit
       ns = ns + 1; ti = tj; ne = ns
       if(ti.gt.rmax) then; ns=ns-1; Exit; end if
      End do

! ... rest of interval with step = hmax

      if(ti.lt.rmax) then
       Do
        ti = ti + hmax
        if(ti.ge.rmax) Exit
        ns = ns + 1
       End do
      End if

      if(allocated(t)) Deallocate(t); Allocate(t(ns+ks))
      t = 0.d0

      Do i=ks+1,ks+ml;       t(i) = hs*(i-ks); End do
      Do i=ks+ml+1,ne;       t(i) = t(i-1) * ht; End do
      if(t(ne).lt.rmax) then
       Do i=ne+1,ns+1;       t(i) = t(i-1) + hmax; End do
      end if
      t(ns+2:ns+ks) = t(ns+1)
      tmax = t(ns+1)

      End Subroutine mkgrid_01


!======================================================================
      Subroutine mkgrid_02
!======================================================================
!     exponential grid for given ns and rmax;  he and hmax are derived.
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_nuclear, z => atomic_number

      Implicit none
      integer :: i
      Real(8) :: ts,ti

      ms=ns+ns; nv=ns-ks+1;  Allocate(t(ns+ks))

! ... initial point:

      if(nuclear.eq.'point') then
       ti=hi/z
      elseif(nuclear.eq.'Fermi') then
       ts = 4.d0 * log(3.d0) * a_fermi
       ti = ts*hi
      else
       ti = r_uniform*hi
      End if

! ... final point:

      tmax=rmax

      he = exp((log(tmax)-log(ti))/(ns-ks))

      t(1:ks)=0.d0; t(ks+1)=ti
      Do i=ks+2,ns+1; t(i)=t(i-1)*he; End do
      t(ns+1:ns+ks)=tmax

      hmax = t(ns+1)-t(ns)

      End Subroutine mkgrid_02



