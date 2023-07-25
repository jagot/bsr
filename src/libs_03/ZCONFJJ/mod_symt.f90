!======================================================================
      Module symt_list
!======================================================================
!     Containes the configuration angular symmetries (i.e. all terms)
!----------------------------------------------------------------------
      Implicit none

      Integer :: nsymt = 0       ! number of symmetries
      Integer :: msymt = 0       ! current dimension of the list
      Integer :: isymt = 10000   ! initial dimension
      Integer :: jsymt = 10      ! average size of one symmetry 
      Integer :: ksymt = 0       ! dimension of all symt.s 
      Integer :: lsymt = 0       ! last element 
      
      Integer,    allocatable :: IT_conf(:)   
      Integer,    allocatable :: ip_term(:)   
      Integer(2), allocatable :: JS_term(:)  
      Integer(1), allocatable :: VS_term(:)  
      Integer(2), allocatable :: JI_term(:)  

! ... IT_stat     - indicate existence the term in given case

      Integer(1), allocatable :: IT_stat(:)

! ... JP_term     - provides ordering of terms according to configurations

      Integer,    allocatable :: JP_term(:)

! ... IT_done(:) - pointer on the done calculation for specific
!                  operators and given terms

      Integer(1), allocatable :: IT_done(:)

! ... IT_need(:) - define the need of calc. for the given term.symmetry

      Integer(1), allocatable :: IT_need(:)

      Integer, parameter ::  mrecl = 100000 ! just to record IT_done

      End Module symt_list


!======================================================================
      Subroutine alloc_symt(m)
!======================================================================
!     allocate, deallocate or reallocate arrays in module symt_list
!----------------------------------------------------------------------
      Use symt_list

      Implicit none
      Integer, Intent(in) :: m
      Integer :: i
      Integer, Allocatable :: iarr(:)

      if(m.le.0) then
       if(allocated(IT_conf)) Deallocate (IT_conf,ip_term,&
                                          JS_term,VS_term,JI_term)
       msymt = 0; nsymt = 0; ksymt = 0; lsymt = 0
      elseif(.not.allocated(IT_conf)) then
       msymt = m; ksymt = msymt*jsymt; lsymt = 0
       Allocate(ip_term(msymt),IT_conf(msymt),&
                JS_term(ksymt),VS_term(ksymt),JI_term(ksymt))
      elseif(m.le.msymt) then
       Return
      elseif(nsymt.eq.0) then
       Deallocate (IT_conf,ip_term,JS_term,VS_term,JI_term)
       msymt = m; ksymt = msymt*jsymt; lsymt = 0
       Allocate(ip_term(msymt),IT_conf(msymt),&
                JS_term(ksymt),VS_term(ksymt),JI_term(ksymt))
      else
       msymt=m; i=lsymt/nsymt+1; if(jsymt.lt.i) jsymt=i; ksymt=jsymt*msymt
       Allocate(iarr(lsymt))
       iarr(1:nsymt)=IT_conf(1:nsymt); Deallocate(IT_conf)
       Allocate(IT_conf(m)); IT_conf=0;  IT_conf(1:nsymt)=iarr(1:nsymt)
       iarr(1:nsymt)=ip_term(1:nsymt); Deallocate(ip_term)
       Allocate(ip_term(m)); ip_term=0;  ip_term(1:nsymt)=iarr(1:nsymt)
       iarr(1:lsymt)=JS_term(1:lsymt); Deallocate(JS_term)
       Allocate(JS_term(ksymt)); JS_term=0;  JS_term(1:lsymt)=iarr(1:lsymt)
       iarr(1:lsymt)=VS_term(1:lsymt); Deallocate(VS_term)
       Allocate(VS_term(ksymt)); VS_term=0;  VS_term(1:lsymt)=iarr(1:lsymt)
       iarr(1:lsymt)=JI_term(1:lsymt); Deallocate(JI_term)
       Allocate(JI_term(ksymt)); JI_term=0;  JI_term(1:lsymt)=iarr(1:lsymt)
       Deallocate(iarr)
      end if

      End Subroutine alloc_symt


!======================================================================
      Integer Function Iadd_symt(iconf,no,Jshell,Vshell,Jintra)
!======================================================================
!     add new overlap conf.symmetry to symt_list
!----------------------------------------------------------------------
      Use symt_list

      Implicit none 
      Integer :: no,iconf, i,j,ip
      Integer, Dimension(*) :: Jshell,Vshell,Jintra

      Iadd_symt = 0
      if(no.le.0) Return
      if(msymt.eq.0) Call Alloc_symt(isymt)

! ... check if the same symt. is already in the list:

      Do i=1,nsymt
       if(iabs(IT_conf(i)).ne.iconf) Cycle
       ip=ip_term(i); Iadd_symt = i
       Do j=1,no; ip = ip+1
        if(Jshell(j).ne.JS_term(ip)) then; Iadd_symt = 0; Exit; end if
        if(Vshell(j).ne.VS_term(ip)) then; Iadd_symt = 0; Exit; end if
        if(Jintra(j).ne.JI_term(ip)) then; Iadd_symt = 0; Exit; end if
       End do
       if(Iadd_symt.ne.0) Return
      End do

! ... Add new symt.:

      if(nsymt.ge.msymt.or.lsymt+no.ge.ksymt) &
         Call Alloc_symt(msymt+isymt)
      nsymt=nsymt+1
      ip_term(nsymt)=lsymt
      IT_conf(nsymt)=iconf
      Do i=1,no; lsymt=lsymt+1
       JS_term(lsymt)=Jshell(i)
       VS_term(lsymt)=Vshell(i)
       JI_term(lsymt)=Jintra(i)
      End do
      Iadd_symt=nsymt

      End Function Iadd_symt


!======================================================================
      Subroutine Get_symt(iterm,iconf,no,Jshell,Vshell,Jintra)
!======================================================================
!     extracts angular symmetry 'iterm'                  
!----------------------------------------------------------------------
      Use symt_list

      Implicit none 
      Integer, Intent(in) :: iterm
      Integer, Intent(out) :: iconf,no
      Integer :: i,ip
      Integer, Dimension(*) :: Jshell,Vshell,Jintra
      Integer, External :: Get_no
      
      if(iterm.le.0.or.iterm.gt.nsymt)  then
        write(*,*) 'Get_symt: it =',iterm, &
                   ' is out of range, nsymt=',nsymt
        Stop
      end if

      iconf=it_conf(iterm)
      no = Get_no(iconf)
      ip = ip_term(iterm)
      Do i=1,no; ip=ip+1
       Jintra(i) = JI_term(ip)
       Jshell(i) = JS_term(ip)
       Vshell(i) = VS_term(ip)
      End do 

      End Subroutine Get_symt


!======================================================================
      Subroutine Read_symt(nu)
!======================================================================
!     reads the list of angular symmetries from unit 'nu'
!----------------------------------------------------------------------
      Use symt_list 

      Implicit none
      Integer :: nu,i

      read(nu) nsymt,lsymt
      if(Allocated(IT_conf)) Deallocate (IT_conf,ip_term, &
                                        JS_term,VS_term,JI_term)
      msymt = nsymt; jsymt=lsymt/nsymt+1; ksymt = lsymt
      Allocate(ip_term(msymt),IT_conf(msymt),&
               JS_term(ksymt),VS_term(ksymt),JI_term(ksymt))
      read(nu) (IT_conf(i),i=1,nsymt)
      read(nu) (ip_term(i),i=1,nsymt)
      read(nu) (JS_term(i),i=1,lsymt)
      read(nu) (VS_term(i),i=1,lsymt)
      read(nu) (JI_term(i),i=1,lsymt)

      End Subroutine Read_symt

!======================================================================
      Subroutine Dummy_symt(nu)
!======================================================================
!     skips the list of angular symmetries in unit 'nu'
!----------------------------------------------------------------------
      Implicit none
      Integer :: nu,i,j,n,l
      Integer(2) :: m
      Integer(1) :: k

      read(nu) n,l
      read(nu) (j,i=1,n)
      read(nu) (j,i=1,n)
      read(nu) (m,i=1,l)
      read(nu) (k,i=1,l)
      read(nu) (m,i=1,l)

      End Subroutine Dummy_symt

!======================================================================
      Subroutine Write_symt(nu)
!======================================================================
!     records the list of angular symmetries to unit 'nu'
!----------------------------------------------------------------------
      Use symt_list 

      Implicit none
      Integer :: nu,i

      write(nu) nsymt,lsymt
      write(nu) (IT_conf(i),i=1,nsymt)
      write(nu) (ip_term(i),i=1,nsymt)
      write(nu) (JS_term(i),i=1,lsymt)
      write(nu) (VS_term(i),i=1,lsymt)
      write(nu) (JI_term(i),i=1,lsymt)

      Deallocate (IT_conf,ip_term,JS_term,VS_term,JI_term)
      msymt = 0; jsymt=lsymt/nsymt+1; ksymt = 0; lsymt=0

      End Subroutine Write_symt


!======================================================================
      Subroutine Write_done(nu)
!======================================================================
!     record the IT_done array to unit "nu"
!---------------------------------------------------------------------- 
      Use symt_list 

      Implicit none
      Integer :: nu,i,i1,i2,n

      n = nsymt*(nsymt+1)/2
      Do i = 1,n
       if(IT_done(i).eq. 0) Stop 'Write_done: IT_done=0 - ?'
       if(IT_done(i).eq.-1) IT_done(i)=0
      End do

      write(nu) n
      i1=1; i2=mrecl; if(i2.gt.n) i2=n
      Do 
       write(nu) (IT_done(i),i=i1,i2)
       i1=i1+mrecl; if(i1.gt.n) Exit
       i2=i2+mrecl; if(i2.gt.n) i2=n
      End do 
 
      End Subroutine Write_done


!======================================================================
      Subroutine Read_done(nu)
!======================================================================
!     reads the IT_done array from unit "nu"
!---------------------------------------------------------------------- 
      Use symt_list 

      Implicit none
      Integer :: nu,i,i1,i2,n

      read(nu) n 
      i1=1; i2=mrecl; if(i2.gt.n) i2=n
      Do 
       read(nu) (IT_done(i),i=i1,i2)
       i1=i1+mrecl; if(i1.gt.n) Exit
       i2=i2+mrecl; if(i2.gt.n) i2=n
      End do 
 
      End Subroutine Read_done


!======================================================================
      Subroutine Dummy_done(nu)
!======================================================================
!     skips  the IT_done array in unit "nu"
!---------------------------------------------------------------------- 
      Use symt_list 

      Implicit none
      Integer :: nu,i,i1,i2,n
      Integer(1) :: j

      read(nu) n 
      i1 = 1; i2 = mrecl; if(i2.gt.n) i2=n
      Do 
       read(nu) (j,i=i1,i2)
       i1=i1+mrecl; if(i1.gt.n) Exit
       i2=i2+mrecl; if(i2.gt.n) i2=n
      End do 
 
      End Subroutine Dummy_done

!======================================================================
      Integer Function Get_iconf(iterm)
!======================================================================
!     provides "iconf" pointer for given "iterm" (for external calls)
!---------------------------------------------------------------------- 
      Use symt_list 
      Implicit none
      Integer, intent(in) :: iterm
      Get_iconf=it_conf(iterm)
      End Function Get_iconf
