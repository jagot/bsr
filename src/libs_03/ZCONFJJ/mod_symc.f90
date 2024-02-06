!======================================================================
      Module symc_list
!======================================================================
!     Containes the "pure-configuration" list  (i.e., without terms)
!----------------------------------------------------------------------
      Implicit none

      Integer :: nsymc = 0       ! number of symmmetries
      Integer :: msymc = 0       ! current dimension of the list
      Integer :: isymc = 10000   ! initial dimension
      Integer :: jsymc = 10      ! average size of one symc.
      Integer :: ksymc = 0       ! max.dimension of all symc.s
      Integer :: lsymc = 0       ! last element

      Integer(1), allocatable :: JT_conf(:)
      Integer(1), allocatable :: no_conf(:)
      Integer,    allocatable :: ip_conf(:)
      Integer(1), allocatable :: iq_conf(:)
      Integer(1), allocatable :: kn_conf(:)

! ... IC_term(ic) - gives for given configuration 'ic' the range of corr.
!                   terms it1,it2 - first and final position in the
!                   ordered list of terms

      Integer, allocatable :: IC_term1(:)
      Integer, allocatable :: IC_term2(:)

! ... IC_need(:) - define the need of calc. between two config.s

      Integer, allocatable :: IC_need(:)

! ... JC_need(:) - define the need of calc. for the given config.

      Integer, allocatable :: JC_need(:)

      End Module symc_list


!======================================================================
      Subroutine alloc_symc(m)
!======================================================================
!     allocate, deallocate, or reallocate arrays in module symc_list
!----------------------------------------------------------------------
      Use symc_list

      Implicit none
      Integer, Intent(in) :: m
      Integer :: i
      Integer, allocatable :: iarr(:)

      if(m.le.0) then
       if(allocated(JT_conf)) Deallocate (JT_conf,no_conf,ip_conf, &
                                          iq_conf,kn_conf)
       msymc = 0; nsymc = 0; ksymc = 0; lsymc = 0
      elseif(.not.allocated(JT_conf)) then
       msymc = m; ksymc = msymc*jsymc; lsymc = 0
       Allocate(JT_conf(msymc),no_conf(msymc),ip_conf(msymc), &
                iq_conf(ksymc),kn_conf(ksymc))
      elseif(m.le.msymc) then
       Return
      elseif(nsymc.eq.0) then
       Deallocate (JT_conf,no_conf,ip_conf,iq_conf,kn_conf)
       msymc = m; ksymc = msymc*jsymc; lsymc = 0
       Allocate(JT_conf(msymc),no_conf(msymc),ip_conf(msymc), &
                iq_conf(ksymc),kn_conf(ksymc))
      else

       msymc=m; i=lsymc/nsymc+1; if(jsymc.lt.i) jsymc=i; ksymc=jsymc*msymc

       Allocate(iarr(lsymc))
       iarr(1:nsymc)=JT_conf(1:nsymc); Deallocate(JT_conf)
       Allocate(JT_conf(m)); JT_conf=0;  JT_conf(1:nsymc)=iarr(1:nsymc)
       iarr(1:nsymc)=no_conf(1:nsymc); Deallocate(no_conf)
       Allocate(no_conf(m)); no_conf=0;  no_conf(1:nsymc)=iarr(1:nsymc)
       iarr(1:nsymc)=ip_conf(1:nsymc); Deallocate(ip_conf)
       Allocate(ip_conf(m)); ip_conf=0;  ip_conf(1:nsymc)=iarr(1:nsymc)
       iarr(1:lsymc)=iq_conf(1:lsymc); Deallocate(iq_conf)
       Allocate(iq_conf(ksymc)); iq_conf=0;  iq_conf(1:lsymc)=iarr(1:lsymc)
       iarr(1:lsymc)=kn_conf(1:lsymc); Deallocate(kn_conf)
       Allocate(kn_conf(ksymc)); kn_conf=0;  kn_conf(1:lsymc)=iarr(1:lsymc)
       Deallocate(iarr)

      end if

      End Subroutine alloc_symc


!======================================================================
      Integer Function Iadd_symc (JT,no,iq,kn)
!======================================================================
!     add new conf.symmetry to the module symc_list
!----------------------------------------------------------------------
      Use symc_list

      Implicit none
      Integer :: no, JT, i,j,ip
      Integer :: iq(*),kn(*)

      Iadd_symc = 0
      if(no.le.0) Return
      if(msymc.le.0) Call Alloc_symc(isymc)

! ... check if the same symc. is already in the list:

      Do i=1,nsymc
       if(JT_conf(i).ne.JT) Cycle
       if(no_conf(i).ne.no) Cycle
       ip=ip_conf(i); Iadd_symc = i
       Do j=1,no; ip=ip+1
        if(iq(j).ne.iq_conf(ip)) then; Iadd_symc = 0; Exit; end if
        if(kn(j).ne.kn_conf(ip)) then; Iadd_symc = 0; Exit; end if
       End do
       if(Iadd_symc.ne.0) Return
      End do

! ... Add new symc.:

      if(nsymc.ge.msymc.or.lsymc+no.ge.ksymc) &
       Call Alloc_symc(msymc+isymc)

      nsymc=nsymc+1
      JT_conf(nsymc)=JT
      no_conf(nsymc)=no
      ip_conf(nsymc)=lsymc
      Do i=1,no; lsymc=lsymc+1
       iq_conf(lsymc)=iq(i)
       kn_conf(lsymc)=kn(i)
      End do
      Iadd_symc=nsymc

      End Function Iadd_symc


!======================================================================
      Subroutine Get_symc(ic,JT,no,nn,kn,ln,jn,iq,in)
!======================================================================
!     extract configuration 'ic' from symc_list
!----------------------------------------------------------------------
      Use symc_list

      Implicit none
      Integer :: ic,JT,no,i,ip
      Integer, Dimension(*) :: nn,kn,ln,jn,iq,in
      Integer, External :: l_kappa, j_kappa

      if(ic.le.0.or.ic.gt.nsymc) Stop 'Get_symc: <ic> is out of range'

      JT = JT_conf(ic)
      no = no_conf(ic)
      ip = ip_conf(ic)
      Do i=1,no; ip = ip+1
       iq(i) = abs(iq_conf(ip))
       kn(i) = kn_conf(ip)
       ln(i) = l_kappa(kn(i))
       jn(i) = j_kappa(kn(i))
       nn(i) = i
       if(iq_conf(ip).lt.0) nn(i)=i-1
       in(i) = 0
      End do

      End Subroutine Get_symc


!======================================================================
      Integer Function Get_no(iconf)
!======================================================================
!     number of shells in configuration 'iconf'
!----------------------------------------------------------------------
      Use symc_list
      Implicit none
      Integer :: iconf
      Get_no = no_conf(iconf)
      End Function Get_no


!======================================================================
      Integer Function Get_jot(iconf)
!======================================================================
!     number of shells in configuration 'iconf'
!----------------------------------------------------------------------
      Use symc_list
      Implicit none
      Integer :: iconf
      Get_jot = JT_conf(iconf)
      End Function Get_jot


!======================================================================
      Subroutine read_symc(nu)
!======================================================================
!     read symc_list from unit "nu"
!----------------------------------------------------------------------
      Use symc_list

      Implicit none
      Integer :: nu,i

      read(nu) nsymc,lsymc
      if(allocated(JT_conf)) Deallocate (JT_conf,no_conf,ip_conf, &
                                         iq_conf,kn_conf)
      Allocate(JT_conf(nsymc),no_conf(nsymc),ip_conf(nsymc), &
               iq_conf(lsymc),kn_conf(lsymc))
      read(nu) (JT_conf(i),i=1,nsymc)
      read(nu) (no_conf(i),i=1,nsymc)
      read(nu) (ip_conf(i),i=1,nsymc)
      read(nu) (iq_conf(i),i=1,lsymc)
      read(nu) (kn_conf(i),i=1,lsymc)
      msymc = nsymc; jsymc = lsymc/nsymc + 1; ksymc = lsymc

      End Subroutine read_symc


!======================================================================
      Subroutine dummy_symc(nu)
!======================================================================
!     dummy read symc_list from unit "nu"
!----------------------------------------------------------------------
      Use symc_list

      Implicit none
      Integer :: nu,i,j,n,l
      Integer(1) :: k

      read(nu) n,l
      read(nu) (k,i=1,n)
      read(nu) (k,i=1,n)
      read(nu) (j,i=1,n)
      read(nu) (k,i=1,l)
      read(nu) (k,i=1,l)

      End Subroutine dummy_symc


!======================================================================
      Subroutine write_symc(nu)
!======================================================================
!     write symc_list to unit "nu"
!----------------------------------------------------------------------
      Use symc_list

      Implicit none
      Integer :: nu,i

      write(nu) nsymc,lsymc
      write(nu) (JT_conf(i),i=1,nsymc)
      write(nu) (no_conf(i),i=1,nsymc)
      write(nu) (ip_conf(i),i=1,nsymc)
      write(nu) (iq_conf(i),i=1,lsymc)
      write(nu) (kn_conf(i),i=1,lsymc)
      Deallocate (JT_conf,no_conf,ip_conf,iq_conf,kn_conf)
      msymc = 0; jsymc = lsymc/nsymc + 1; ksymc = 0; lsymc=0

      End Subroutine write_symc
