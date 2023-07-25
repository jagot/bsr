!======================================================================
      Module symt_list_LS
!======================================================================
!
!     Containes the configuration symmetries.
!
!     ROUTINES:  alloc_symt_LS  (m)
!                Iadd_symt_LS   (iconf,no,LS)
!                Get_symt_LS    (it,iconf,no,LS)
!                Write_symt_LS  (nu)
!                Read_symt_LS   (nu)
!                Write_oper_LS  (nu)
!                Read_oper_LS   (nu)
!                no_term_LS     (ic)      
!----------------------------------------------------------------------
      Use param_LS

      Implicit none

      Integer :: nsymt = 0       ! number of symt.s
      Integer :: msymt = 0       ! current dimension of the list
      Integer :: isymt = 2**12   ! initial dimension
      Integer :: jsymt = 2*3     ! average size of one symt. 
      Integer :: ksymt = 0       ! dimension of all symt.s 
      Integer :: lsymt = 0       ! last element 
      
      Integer, allocatable :: IT_conf(:)   ! configuration pointer
      Integer, allocatable :: ip_term(:)   ! position pointer
      Integer, allocatable :: LS_term(:,:) ! main array of LS terms  

! ... IT_stat     - additonal index

      Integer, allocatable :: IT_stat(:)

! ... IT_sort     - provides ordering of terms according to configurations

      Integer, allocatable :: IT_sort(:)

! ... IT_done(:) - pointer on the done calculation for specific
!                  operators and given terms

      Integer(1), allocatable :: IT_done(:)

      Integer(1), allocatable :: IT_oper(:,:)


! ... IT_need(:) - define the need of calc. for the given term.symmetry

      Integer(1), allocatable :: IT_need(:), JT_need(:)

      Real(8) :: mem_symt = 0.d0,  mem_oper = 0.d0
      Integer :: m_symt

      Integer(8) :: ij, ij_oper, ii8, jj8

      End Module symt_list_LS


!======================================================================
      Subroutine alloc_symt_LS(m)
!======================================================================
!     allocate arrays in the module "symt_list_LS"
!----------------------------------------------------------------------
      Use symt_list_LS

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: ia(:),ib(:,:)

      if(m.le.0) then
       if(allocated(IT_conf)) Deallocate (IT_conf,LS_term,ip_term,IT_stat)
       msymt = 0; nsymt = 0; ksymt = 0; lsymt = 0
      elseif(.not.allocated(IT_conf)) then
       msymt = m; ksymt = msymt*jsymt; lsymt = 0
       Allocate(IT_conf(msymt),ip_term(msymt),LS_term(ksymt,5),IT_stat(msymt))
       IT_stat=0; IT_conf=0; ip_term=0; LS_term=0
      elseif(m.le.msymt) then
       Return
      elseif(nsymt.eq.0) then
       Deallocate (IT_conf,ip_term,LS_term,IT_stat)
       msymt = m; ksymt = msymt*jsymt; lsymt = 0
       Allocate(IT_conf(msymt),ip_term(msymt),LS_term(ksymt,5),IT_stat(msymt))
       IT_stat=0; IT_conf=0; ip_term=0; LS_term=0
      else
       msymt = m; jsymt=lsymt/nsymt+1; ksymt = msymt*jsymt
       if(ksymt.lt.lsymt) ksymt=lsymt
	Allocate(ia(nsymt))
	ia=IT_conf(1:nsymt); Deallocate(IT_conf)
	Allocate(IT_conf(msymt)); IT_conf(1:nsymt)=ia
	ia=ip_term(1:nsymt); Deallocate(ip_term)
	Allocate(ip_term(msymt)); ip_term(1:nsymt)=ia
	ia=IT_stat(1:nsymt); Deallocate(IT_stat)
	Allocate(IT_stat(msymt)); IT_stat=0; IT_stat(1:nsymt)=ia
       Deallocate(ia)
	Allocate(ib(lsymt,5))
	ib=LS_term(1:lsymt,:); Deallocate(LS_term)
	Allocate(LS_term(ksymt,5)); LS_term(1:lsymt,:)=ib
       Deallocate(ib)
      end if

      mem_symt = (3.*msymt + 5.*ksymt)*4/(1024*1024) 
      m_symt = 3*msymt + 5*ksymt

      End Subroutine alloc_symt_LS


!======================================================================
      Subroutine alloc_it_oper_LS(m)
!======================================================================
      Use symt_list_LS

      Implicit none
      Integer, intent(in) :: m

      if(allocated(it_oper)) Deallocate(it_oper)
      mem_oper=0.d0       
      if(m.le.0) Return

      ij = nsymt; ij_oper =  ij*(ij+1)/2
      Allocate(it_oper(noper,ij_oper))
      it_oper = 0
      mem_oper = 1.d0*(noper+1)*ij_oper/(1024*1024)

      End Subroutine alloc_it_oper_LS


!======================================================================
      Integer Function Iadd_symt_LS(iconf,no,LS)
!======================================================================
!     add new overlap conf.symmetry to symt_list_LS
!----------------------------------------------------------------------
      Use symt_list_LS

      Implicit none 
      Integer :: iconf,no, i,j,ip,jp
      Integer :: LS(msh,5)
      Integer, external :: no_conf_LS

      Iadd_symt_LS = 0
      if(no.le.0) Return
      if(msymt.eq.0) Call Alloc_symt_LS(isymt)
      if(no.ne.no_conf_LS(iconf)) Stop 'Iadd_symt_LS: no <> no_conf'

! ... check if the same symt. is already in the list:

      Do i=1,nsymt
       if(iabs(IT_conf(i)).ne.iconf) Cycle
       ip=ip_term(i); Iadd_symt_LS = i
       Do j=1,no; ip=ip+1
        if(LS(j,1).ne.LS_term(ip,1)) then; Iadd_symt_LS=0; Exit; end if
        if(LS(j,2).ne.LS_term(ip,2)) then; Iadd_symt_LS=0; Exit; end if
        if(LS(j,3).ne.LS_term(ip,3)) then; Iadd_symt_LS=0; Exit; end if
        if(LS(j,4).ne.LS_term(ip,4)) then; Iadd_symt_LS=0; Exit; end if
        if(LS(j,5).ne.LS_term(ip,5)) then; Iadd_symt_LS=0; Exit; end if
       End do
       if(Iadd_symt_LS.ne.0) Return
      End do

! ... Add new symt:

      if(nsymt.ge.msymt.or.lsymt+no.gt.ksymt) &
      Call Alloc_symt_LS(msymt+isymt)
      nsymt=nsymt+1
      IT_conf(nsymt)=iconf
      ip_term(nsymt)=lsymt
      Do i=1,no; lsymt=lsymt+1; LS_term(lsymt,:)=LS(i,:); End do

      Iadd_symt_LS=nsymt

      End Function Iadd_symt_LS


!======================================================================
      Subroutine Get_symt_LS(it,iconf,no,LS)
!======================================================================
!     extrats term symetry 'it'                  
!----------------------------------------------------------------------
      Use symt_list_LS

      Implicit none 
      Integer :: it,iconf, no, i,j,ip
      Integer :: LS(msh,5)
      Integer, external :: no_conf_LS
      
      if(it.le.0.or.it.gt.nsymt) Stop 'Get_symt_LS: <it> is out of range'

      iconf = IT_conf(it)
      no = no_conf_LS(iconf)
      ip = ip_term(it)
      Do i=1,no; ip=ip+1;  LS(i,:) = LS_term(ip,:);  End do 

      End Subroutine Get_symt_LS


!======================================================================
      Subroutine Read_symt_LS(nu)
!======================================================================
      Use symt_list_LS 

      Implicit none
      Integer :: nu,i

      read(nu) nsymt,lsymt
      if(Allocated(IT_conf)) Deallocate (IT_conf,ip_term,LS_term,IT_stat)
      msymt = nsymt; jsymt=lsymt/nsymt+1; ksymt = lsymt
      Allocate(ip_term(msymt),IT_conf(msymt),LS_term(ksymt,5),IT_stat(msymt))
      read(nu) (IT_conf(i),i=1,nsymt)
      read(nu) (ip_term(i),i=1,nsymt)
      read(nu) (LS_term(i,:),i=1,lsymt)
      IT_stat = 0

      mem_symt = (3.*msymt + 5.*ksymt)*4/(1024*1024) 

      End Subroutine Read_symt_LS


!======================================================================
      Subroutine Write_symt_LS(nu)
!======================================================================
      Use symt_list_LS 

      Implicit none
      Integer :: nu,i

      write(nu) nsymt,lsymt
      write(nu) (IT_conf(i),i=1,nsymt)
      write(nu) (ip_term(i),i=1,nsymt)
      write(nu) (LS_term(i,:),i=1,lsymt)
      Deallocate (IT_conf,ip_term,LS_term,IT_stat)
      msymt = 0; jsymt=lsymt/nsymt+1; ksymt = 0; lsymt=0

      End Subroutine Write_symt_LS

!======================================================================
      Subroutine Write_oper_LS(nu)
!======================================================================
      Use symt_list_LS 

      Implicit none
      Integer :: nu,i,j,n

      n = nsymt*(nsymt+1)/2
      write(nu) n
      Do i = 1,n
       Do j=1,noper
        if(IT_oper(j,i).eq. 0) IT_oper(j,i)=1
        if(IT_oper(j,i).eq.-1) IT_oper(j,i)=0
       End do 
       write(nu) IT_oper(:,i)
      End do

      End Subroutine Write_oper_LS


!======================================================================
      Subroutine Record_oper_LS(nu)
!======================================================================
      Use symt_list_LS 

      Implicit none
      Integer, intent(in) :: nu
      Integer :: j

      write(nu) ij_oper
      Do ij = 1,ij_oper
       Do j=1,noper
        if(IT_oper(j,ij).eq. 0) IT_oper(j,ij)=1
        if(IT_oper(j,ij).eq.-1) IT_oper(j,ij)=0
       End do 
      End do
      write(nu) IT_oper

      End Subroutine Record_oper_LS


!======================================================================
      Subroutine Load_oper_LS(nu)
!======================================================================
      Use symt_list_LS 

      Implicit none
      Integer, intent(in) :: nu
      Integer(8) :: i,j

      read(nu) ij
      if(ij_oper.lt.ij) Stop 'Read_oper_LS: nsymt too small' 
      read(nu) ((IT_oper(i,j),i=1,noper),j=1,ij)

      End Subroutine Load_oper_LS


!======================================================================
      Subroutine Read_oper_LS(nu)
!======================================================================
      Use symt_list_LS 

      Implicit none
      Integer :: nu,i,n

      read(nu) n
      if(nsymt*(nsymt+1)/2.lt.n) Stop 'Read_oper_LS: nsymt too small' 
      if(.not.Allocated(IT_oper)) Allocate(IT_oper(noper,n))
      IT_oper=0
      Do i = 1,n; read(nu) IT_oper(:,i); End do

      End Subroutine Read_oper_LS


!======================================================================
      Subroutine Write_done_LS(nu)
!======================================================================
      Use symt_list_LS, only: nsymt, IT_done 

      Implicit none
      Integer :: nu,i,i1,i2,n
      Integer, parameter :: m=1000

      n = nsymt*(nsymt+1)/2
      write(nu) n
      i1=1; i2=m; if(i2.gt.n) i2=n
    1 write(nu) (IT_done(i),i=i1,i2)
      i1=i1+m; i2=i2+m; if(i2.gt.n) i2=n
      if(i1.le.n) go to 1

      End Subroutine Write_done_LS

!======================================================================
      Subroutine Read_done_LS(nu)
!======================================================================
      Use symt_list_LS, only: IT_done 

      Implicit none
      Integer :: nu,i,i1,i2,n
      Integer, parameter :: m=1000

      read(nu) n
      if(.not.Allocated(IT_done)) then
        Allocate(IT_done(n))
      else
        i=SIZE(IT_done)
        if(i.lt.n) then
         Deallocate(IT_done); Allocate(IT_done(n))
        end if 
      end if

      i1=1; i2=m; if(i2.gt.n) i2=n
    1 read(nu) (IT_done(i),i=i1,i2)
      i1=i1+m; i2=i2+m; if(i2.gt.n) i2=n
      if(i1.le.n) go to 1

      End Subroutine Read_done_LS


!======================================================================
      Integer Function no_term_LS (it)
!======================================================================
!     number of shells in state 'iconf'                   
!----------------------------------------------------------------------
      Use symt_list_LS

      Implicit none 
      Integer :: it
      Integer, External :: no_conf_LS

      no_term_LS = no_conf_LS(it_conf(it))

      End Function no_term_LS

