!======================================================================
      Module symc_list_LS
!======================================================================
!     Containes the configuration symmetries, defined as description
!     of configurations without principal quantum numbers
!
!     ROUTINES:   alloc_symc_LS (m)
!                 Iadd_symc_LS  (LT,ST,no,iq,ln)
!                 Get_symc_LS   (ic,LT,ST,no,nn,ln,iq,kn)
!                 read_symc_LS  (nu)
!                 write_symc_LS (nu)
!                 def_maxl      (l)
!                 no_conf_LS    (ic)
!----------------------------------------------------------------------
      Implicit none

      Integer :: nsymc = 0       ! number of symmmetries
      Integer :: msymc = 0       ! current dimension of the list
      Integer :: isymc = 2**12   ! initial dimension
      Integer :: jsymc = 2*3     ! average size of one symc.
      Integer :: ksymc = 0       ! max.dimension of all symc.s
      Integer :: lsymc = 0       ! last element

      Integer, allocatable :: LT_conf(:), ST_conf(:)  ! total terms
      Integer, allocatable :: no_conf(:)              ! number of shells
                                                      ! (we can avoid this array)
      Integer, allocatable :: ip_conf(:)              ! pointer for the given config.
                                                      ! in the common list {iq,ln}
      Integer, allocatable :: iq_conf(:)              ! occupation numbers
      Integer, allocatable :: ln_conf(:)              ! l-values

! ... IC_term1(ic), IC_term2(ic) - gives, for the configuration 'ic',
! ... the range of related terms in the ordered list of terms

      Integer, allocatable :: IC_term1(:), IC_term2(:)

! ... IC_need(:) - define the need of calc. for two given config.s
! ... JC_need(:) - define the need of calc. for the given config.

      Integer, allocatable :: IC_need(:), JC_need(:)

      Integer :: m_symc ! memory requirements

      End Module symc_list_LS


!======================================================================
      Subroutine alloc_symc_LS(m)
!======================================================================
!     allocate arrays in the module "symc_list_LS"
!----------------------------------------------------------------------
      Use symc_list_LS

      Implicit none
      Integer, Intent(in) :: m
      Integer, allocatable :: ia(:)

      if(m.le.0) then
       if(allocated(LT_conf)) Deallocate (LT_conf,no_conf,ip_conf, &
                                          iq_conf,ST_conf,ln_conf)
       msymc = 0; nsymc = 0; ksymc = 0; lsymc = 0
      elseif(.not.allocated(LT_conf)) then
       msymc = m; ksymc = msymc*jsymc; lsymc = 0
       Allocate(LT_conf(msymc),ST_conf(msymc),no_conf(msymc), &
                ip_conf(msymc),iq_conf(ksymc),ln_conf(ksymc))
      elseif(m.le.msymc) then
       Return
      elseif(nsymc.eq.0) then
       Deallocate (LT_conf,no_conf,ip_conf,iq_conf,ST_conf,ln_conf)
       msymc = m; ksymc = msymc*jsymc; lsymc = 0
       Allocate(LT_conf(msymc),ST_conf(msymc),no_conf(msymc), &
                ip_conf(msymc),iq_conf(ksymc),ln_conf(ksymc))
      else
       msymc = m; jsymc = lsymc/nsymc + 1; ksymc = jsymc*msymc
       if(ksymc.lt.lsymc) ksymc=lsymc
       Allocate(ia(nsymc))
       ia = LT_conf(1:nsymc); Deallocate(LT_conf)
       Allocate(LT_conf(msymc)); LT_conf(1:nsymc)=ia
       ia = ST_conf(1:nsymc); Deallocate(ST_conf)
       Allocate(ST_conf(msymc)); ST_conf(1:nsymc)=ia
       ia = no_conf(1:nsymc); Deallocate(no_conf)
       Allocate(no_conf(msymc)); no_conf(1:nsymc)=ia
       ia = ip_conf(1:nsymc); Deallocate(ip_conf)
       Allocate(ip_conf(msymc)); ip_conf(1:nsymc)=ia
       Deallocate(ia)
       Allocate(ia(lsymc))
       ia = iq_conf(1:lsymc); Deallocate(iq_conf)
       Allocate(iq_conf(ksymc)); iq_conf(1:lsymc)=ia
       ia = ln_conf(1:lsymc); Deallocate(ln_conf)
       Allocate(ln_conf(ksymc)); ln_conf(1:lsymc)=ia
       Deallocate(ia)
      end if

      m_symc = 4*msymc + 2*ksymc

      End Subroutine alloc_symc_LS


!======================================================================
      Integer Function Iadd_symc_LS (LT,ST,no,iq,ln)
!======================================================================
!     add new overlap conf.symmetry to symc_list
!----------------------------------------------------------------------
      Use symc_list_LS

      Implicit none
      Integer :: no, LT,ST,  i,j,ip
      Integer, Dimension(no) :: iq,ln

      Iadd_symc_LS = 0
      if(no.le.0) Return
      if(msymc.eq.0) Call Alloc_symc_LS(isymc)

! ... check if the same symc. is already in the list:

      Do i=1,nsymc
       if(LT_conf(i).ne.LT) Cycle
       if(ST_conf(i).ne.ST) Cycle
       if(no_conf(i).ne.no) Cycle
       ip=ip_conf(i); Iadd_symc_LS = i
       Do j=1,no; ip=ip+1
        if(iq(j).ne.iq_conf(ip)) then; Iadd_symc_LS=0; Exit; end if
        if(ln(j).ne.ln_conf(ip)) then; Iadd_symc_LS=0; Exit; end if
       End do
       if(Iadd_symc_LS.ne.0) Return
      End do

! ... Add new symc.:

      if(nsymc.ge.msymc.or.lsymc+no.gt.ksymc) &
       Call Alloc_symc_LS(msymc+isymc)

      nsymc=nsymc+1
      LT_conf(nsymc)=LT
      ST_conf(nsymc)=ST
      no_conf(nsymc)=no
      ip_conf(nsymc)=lsymc
      Do i=1,no; lsymc=lsymc+1
       iq_conf(lsymc)=iq(i)
       ln_conf(lsymc)=ln(i)
      End do
      Iadd_symc_LS = nsymc

      End Function Iadd_symc_LS


!======================================================================
      Subroutine Get_symc_LS(ic,LT,ST,no,nn,ln,iq,kn)
!======================================================================
!     extract configuration 'ic'
!----------------------------------------------------------------------
      Use symc_list_LS

      Implicit none
      Integer :: ic,LT,ST,no,i,ip
      Integer :: nn(*),ln(*),iq(*),kn(*)

      if(ic.le.0.or.ic.gt.nsymc) Stop 'Get_symc: <ic> is out of range'

      LT = LT_conf(ic)
      ST = ST_conf(ic)
      no = no_conf(ic)
      ip = ip_conf(ic)
      Do i=1,no; ip=ip+1
       iq(i) = iq_conf(ip)
       ln(i) = ln_conf(ip)
       nn(i) = i
       kn(i) = 0
      End do

      End Subroutine Get_symc_LS


!======================================================================
      Subroutine read_symc_LS(nu)
!======================================================================
      Use symc_list_LS

      Implicit none
      Integer :: nu,i

      read(nu) nsymc,lsymc
      if(allocated(LT_conf)) &
      Deallocate (LT_conf,ST_conf,no_conf,ip_conf,iq_conf,ln_conf)
      Allocate(LT_conf(nsymc),ST_conf(lsymc),no_conf(nsymc),ip_conf(nsymc), &
               iq_conf(lsymc),ln_conf(lsymc))
      read(nu) (LT_conf(i),i=1,nsymc)
      read(nu) (ST_conf(i),i=1,nsymc)
      read(nu) (no_conf(i),i=1,nsymc)
      read(nu) (ip_conf(i),i=1,nsymc)
      read(nu) (iq_conf(i),i=1,lsymc)
      read(nu) (ln_conf(i),i=1,lsymc)
      msymc = nsymc; jsymc = lsymc/nsymc + 1; ksymc = lsymc

      End Subroutine read_symc_LS


!======================================================================
      Subroutine write_symc_LS (nu)
!======================================================================
      Use symc_list_LS

      Implicit none
      Integer :: nu,i

      write(nu) nsymc,lsymc
      write(nu) (LT_conf(i),i=1,nsymc)
      write(nu) (ST_conf(i),i=1,nsymc)
      write(nu) (no_conf(i),i=1,nsymc)
      write(nu) (ip_conf(i),i=1,nsymc)
      write(nu) (iq_conf(i),i=1,lsymc)
      write(nu) (ln_conf(i),i=1,lsymc)
      Deallocate (LT_conf,ST_conf,no_conf,ip_conf,iq_conf,ln_conf)
      msymc = 0; jsymc = lsymc/nsymc + 1; ksymc = 0; lsymc=0

      End Subroutine write_symc_LS


!======================================================================
      Subroutine Def_maxl(l)
!======================================================================
      Use symc_list_LS

      Implicit none
      Integer :: l,i

      l=0;  Do i=1,lsymc; if(l.lt.ln_conf(i)) l=ln_conf(i); End do

      End Subroutine Def_maxl


!======================================================================
      Integer Function no_conf_LS (ic)
!======================================================================
!     number of shells in state 'iconf'
!----------------------------------------------------------------------
      Use symc_list_LS

      Implicit none
      Integer :: ic

      no_conf_LS = no_conf(ic)

      End Function no_conf_LS



