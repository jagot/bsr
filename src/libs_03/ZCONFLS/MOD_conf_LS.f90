!=====================================================================
      Module conf_LS
!=====================================================================
!     description of configuration lists
!     ROUTINES:  alloc_cfg_LS (m)
!                Iadd_cfg_LS  ()
!                Ifind_cfg_LS ()
!                Jfind_cfg_LS ()
!                Get_cfg_LS   (ic)
!                Save_cfg_LS  (i)
!                no_ic_LS     (ic)
!---------------------------------------------------------------------
      Use param_LS

      Implicit none
      Integer :: ne     = 0     !  number of electrons
      Integer :: parity = 0     !  parity of states (+1,-1)

! ... description of 1 conf.w.function:

      Integer :: no, Ltotal, Stotal, Jtotal, Ptotal, iconf, iterm
      Integer, dimension(msh) :: nn,ln,iq,kn, np_orb, np
      Integer :: LS(msh,5)

      Integer :: no1, Ltotal1,Stotal1,Jtotal1,Ptotal1, iconf1,iterm1
      Integer, dimension(msh) :: nn1,ln1,iq1,kn1, np_orb1, np1
      Integer :: LS1(msh,5)

      Integer :: no2, Ltotal2,Stotal2,Jtotal2,Ptotal2, iconf2,iterm2
      Integer, dimension(msh) :: nn2,ln2,iq2,kn2, np_orb2, np2
      Integer :: LS2(msh,5)

      Integer, allocatable :: LSI(:,:,:)

! ... Storing configurations as character strings:

      Character(200) :: CONFIG, COUPLE

! ... closed shells core

      Real(8) :: Ecore = 0.d0
      Integer :: nclosd = 0, ncore = 0
      Character(200) :: closed='    ', core='    '

! ... CONFIGURATION LIST PARAMETERS

      Integer :: ncfg  = 0       !  current number of configurations
      Integer :: mcfg  = 0       !  max. number of configurations
      Integer :: icfg  = 1000000 !  initial prediction of mcfg
      Integer :: jcfg  = 5       !  average number of shells
      Integer :: kcfg  = 0       !  max. dimensionr
      Integer :: lcfg  = 0       !  last element

      Integer :: nzero = 0       !  first-order set in config. list
      Integer :: ncfg1,ncfg2
      Integer :: kshift1,kshift2

      Integer :: L_min,L_max, S_min,S_max, J_min,J_max, conf_mode = 0

! ... expansion coeficients

      Real(8), allocatable :: WC(:)

! ... label representation of configurations

      Character(200), allocatable :: LABEL(:)

! ... configuration overlap tolerance:

      Real(8) :: c_ovl = -1.d0

! ... additional arrays:

      Integer, allocatable :: IP_state(:)
      Integer, allocatable :: IC_term(:)
      Integer, allocatable :: IP_orb(:)

      Integer, allocatable :: IP_stat(:)
      Integer, allocatable :: IT_state1(:)
      Integer, allocatable :: IT_state2(:)
      Integer, allocatable :: IC_sym(:)

      Integer :: m_conf_LS

      End Module conf_LS


!======================================================================
      Subroutine alloc_cfg_LS(m)
!======================================================================
! ... allocate the configurations in module "conf_LS"
!----------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer :: m
      Integer, allocatable :: ia(:)
      Real(8), allocatable :: ra(:)

      if(m.le.0) then
       if(allocated(ip_state)) Deallocate (ip_state,ip_orb,IC_term,WC)
       mcfg = 0; ncfg = 0; lcfg = 0; kcfg = 0
      elseif(.not.allocated(ip_state)) then
       mcfg = m; kcfg = mcfg*jcfg
       Allocate(ip_state(mcfg),IC_term(mcfg),ip_orb(kcfg),WC(mcfg))
       ncfg = 0; lcfg = 0; WC = 0.d0
      elseif(m.le.mcfg.and.lcfg+no.lt.kcfg) then
       Return
      elseif(ncfg.le.0) then
       Deallocate (ip_state,ip_orb,IC_term,WC); ncfg = 0; lcfg = 0
       mcfg = m; kcfg = mcfg*jcfg
       Allocate(ip_state(mcfg),IC_term(mcfg),ip_orb(kcfg),WC(mcfg))
       WC = 0.d0
      else
       mcfg=m; jcfg = lcfg/ncfg + 1; kcfg=mcfg*jcfg
       if(kcfg.lt.lcfg) kcfg=lcfg
       Allocate(ia(ncfg))
       ia=IC_term(1:ncfg); Deallocate(IC_term)
       Allocate(IC_term(mcfg)); IC_term(1:ncfg)=ia
       ia=ip_state(1:ncfg); Deallocate(ip_state)
       Allocate(ip_state(mcfg)); ip_state(1:ncfg)=ia
       Deallocate(ia)
       Allocate(ia(lcfg))
       ia=ip_orb(1:lcfg); Deallocate(ip_orb)
       Allocate(ip_orb(kcfg)); ip_orb(1:lcfg)=ia
       Deallocate(ia)
       Allocate(ra(ncfg))
       ra=WC(1:ncfg); Deallocate(WC); Allocate(WC(mcfg));
       WC=0.d0; WC(1:ncfg)=ra
       Deallocate(ra)
      end if

      m_conf_LS = 4*mcfg + kcfg +  200 + 4*mcfg   ! ??? last is optional

      End Subroutine alloc_cfg_LS


!======================================================================
      Integer Function Iadd_cfg_LS()
!======================================================================
!     add new CAS to cfg_list
!----------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer :: i,ic,ip
      Integer, external :: Iadd_symc_LS, Iadd_symt_LS, Ifind_nlk

      Iadd_cfg_LS = 0
      if(no.le.0) Return
      if(mcfg.eq.0) Call Alloc_cfg_LS(icfg)

      iconf = Iadd_symc_LS(LS(no,4),LS(no,5),no,iq,ln)
      iterm = Iadd_symt_LS(iconf,no,LS)
      Do i = 1,no; np(i)=Ifind_nlk(nn(i),ln(i),kn(i),2); End do

      ncfg=ncfg+1
      if(ncfg.eq.mcfg.or.lcfg+no.gt.kcfg) Call Alloc_cfg_LS(mcfg+icfg)

      IC_term(ncfg)=iterm
      ip_state (ncfg)=lcfg
      Do i=1,no; lcfg=lcfg+1; IP_orb(lcfg)=np(i); End do
      Iadd_cfg_LS = ncfg

      End Function Iadd_cfg_LS


!======================================================================
      Integer Function Ifind_cfg_LS()
!======================================================================
!     find or add new CAS to cfg_list
!----------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer :: i,ic,ip
      Integer, External :: Iadd_symc_LS, Iadd_symt_LS, Ifind_nlk

      Ifind_cfg_LS = 0
      if(no.le.0) Return

      if(mcfg.eq.0) Call Alloc_cfg_LS(icfg)

      iconf = Iadd_symc_LS(LS(no,4),LS(no,5),no,iq,ln)
      iterm = Iadd_symt_LS(iconf,no,LS)
      Do i = 1,no; np(i)=Ifind_nlk(nn(i),ln(i),kn(i),2); End do

! ... check if we already have such state:

      Do ic = 1,ncfg
       if(IC_term(ic).ne.iterm) Cycle
       ip = ip_state(ic); Ifind_cfg_LS = ic
       Do i = 1,no; ip=ip+1
        if(np(i).ne.IP_orb(ip)) then; Ifind_cfg_LS=0; Exit; end if
       End do
       if(Ifind_cfg_LS.ne.0) Return
      End do

      ncfg=ncfg+1
      if(ncfg.eq.mcfg.or.lcfg+no.gt.kcfg) Call Alloc_cfg_LS(mcfg+icfg)

      IC_term(ncfg)=iterm
      ip_state (ncfg)=lcfg
      Do i=1,no; lcfg=lcfg+1; IP_orb(lcfg)=np(i); End do
      Ifind_cfg_LS = ncfg

      End Function Ifind_cfg_LS


!======================================================================
      Integer Function Jfind_cfg_LS()
!======================================================================
!     find or add new CAS to cfg_list
!     Jfind_cfg_LS = configuration index (<0 if was found)
!----------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer :: i,ic,ip
      Integer, External :: Iadd_symc_LS, Iadd_symt_LS, Ifind_nlk

      Jfind_cfg_LS = 0
      if(no.le.0) Return

      if(mcfg.eq.0) Call Alloc_cfg_LS(icfg)

      iconf = Iadd_symc_LS(LS(no,4),LS(no,5),no,iq,ln)
      iterm = Iadd_symt_LS(iconf,no,LS)
      Do i = 1,no; np(i)=Ifind_nlk(nn(i),ln(i),kn(i),2); End do

! ... check if we already have such state:

      Do ic = 1,ncfg
       if(IC_term(ic).ne.iterm) Cycle
       ip = ip_state(ic); Jfind_cfg_LS = -ic
       Do i = 1,no; ip=ip+1
        if(np(i).ne.IP_orb(ip)) then; Jfind_cfg_LS=0; Exit; end if
       End do
       if(Jfind_cfg_LS.ne.0) Return
      End do

      ncfg=ncfg+1
      if(ncfg.eq.mcfg.or.lcfg+no.gt.kcfg) Call Alloc_cfg_LS(mcfg+icfg)

      IC_term(ncfg)=iterm
      ip_state (ncfg)=lcfg
      Do i=1,no; lcfg=lcfg+1; IP_orb(lcfg)=np(i); End do
      Jfind_cfg_LS = ncfg

      End Function Jfind_cfg_LS


!======================================================================
      Subroutine Get_cfg_LS(ic)
!======================================================================
!     extract the configuration (ic) from the cfg-list
!----------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer :: ic, i,j,ip, l,ii

      no=0; nn=0; ln=0; kn=0; LS=0
      iterm=IC_term(ic)

      Call Get_symt_LS(iterm,iconf,no,LS)

      Call Get_symc_LS(iconf,Ltotal,Stotal,no,nn,ln,iq,kn)

      if(Ltotal.ne.LS(no,4).or.Stotal.ne.LS(no,5)) then
       write(*,*) 'Total term in symc is not consistent with symt:'
       write(*,*) 'ic, iconf, iterm =',ic,iconf,iterm
       write(*,*) 'Ltotal,Stotal,LS(no) = ', Ltotal,Stotal,LS(no,4:5)
       Stop ' '
      end if
      ip = ip_state(ic)
      Do i=1,no; ip=ip+1
       j = IP_orb(ip); Call Get_nlki (j,nn(i),l,kn(i),ii)
      End do

      End Subroutine Get_cfg_LS


!======================================================================
      Subroutine Save_cfg_LS(i)
!======================================================================
!     save(restore) curent state in position i
!----------------------------------------------------------------------
      Use conf_LS
      Implicit none
      Integer, intent(in) :: i

      Select Case(i)
      Case(1)
       no1=no; nn1=nn; kn1=kn; ln1=ln; iq1=iq
       LS1 = LS; Stotal1=Stotal; Ltotal1=Ltotal
       iconf1=iconf; iterm1=iterm
      Case(2)
       no2=no; nn2=nn; kn2=kn; ln2=ln; iq2=iq
       LS2 = LS; Stotal2=Stotal; Ltotal2=Ltotal
       iconf2=iconf; iterm2=iterm
      Case(-1)
       no=no1; nn=nn1; kn=kn1; ln=ln1; iq=iq1
       LS=LS1; Stotal=Stotal1; Ltotal=Ltotal1
       iconf=iconf1; iterm=iterm1
      Case(-2)
       no=no2; nn=nn2; kn=kn2; ln=ln2; iq=iq2
       LS=LS2; Stotal=Stotal2; Ltotal=Ltotal2
       iconf=iconf2; iterm=iterm2
       Case Default
       Stop 'save_cfg: ???'
      End Select

      End Subroutine Save_cfg_LS


!======================================================================
      Integer Function no_ic_LS (ic)
!======================================================================
!     number of shells in state 'ic'
!----------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer :: ic
      Integer, External :: no_term_LS

      no_ic_LS = no_term_LS ( IC_term(ic) )

      End Function no_ic_LS


!======================================================================
      Subroutine Write_conf_LS(nu)
!======================================================================
      Use conf_LS

      Implicit none
      Integer :: nu,i

      write(nu) ncfg,lcfg,ne
      write(nu) (ip_state(i),i=1,ncfg)
      write(nu) (ic_term(i),i=1,ncfg)
      write(nu) (ip_orb(i),i=1,lcfg)
      write(nu) (WC(i),i=1,ncfg)

      End Subroutine Write_conf_LS


!======================================================================
      Subroutine Read_conf_LS(nu)
!======================================================================
      Use conf_LS

      Implicit none
      Integer :: nu,i

      Call Alloc_cfg_LS(0)
      read(nu) ncfg,lcfg,ne
      Allocate(ip_state(ncfg),IC_term(ncfg),ip_orb(lcfg),WC(ncfg))
      mcfg = ncfg; kcfg = lcfg; jcfg = kcfg/mcfg + 1
      read(nu) (ip_state(i),i=1,ncfg)
      read(nu) (ic_term(i),i=1,ncfg)
      read(nu) (ip_orb(i),i=1,lcfg)
      read(nu) (WC(i),i=1,ncfg)

      End Subroutine Read_conf_LS


!======================================================================
      Subroutine Symc_conf(ic,conf)
!======================================================================
!     get label for configuration
!----------------------------------------------------------------------
      Use conf_LS, only: LTOTAL,STOTAL,no,nn,ln,iq,kn

      Implicit none
      Integer, intent(in) :: ic
      Integer :: i,k
      Character(*) :: conf
      Character(4), external :: AL

      conf = ' '

      Call Get_symc_LS(ic,LTOTAL,STOTAL,no,nn,ln,iq,kn)

      k=0
      Do i=1,no; if(iq(i).eq.0) Cycle
       conf(k+2:k+2) = AL(ln(i),1)
       write(conf(k+3:k+6),'(a1,i2,a1)') '(',iq(i),')'
       k=k+6
      End do

      End Subroutine Symc_conf


