!=====================================================================
      Module conf_jj
!=====================================================================
!     containes description of the configuration atomic states
!     and related information
!
!     The configuration atomic states (CAS) are defined by "iterm"
!     pointer on the "angular symmetry" and by list of specific
!     orbitals, which are placed in array ip_orb:
!     iterm     => IS_term(ic);
!     orbitals  => ip_orb(ip+1,...ip+no), where ip => IP_state(ic)
!     ip_orb(i) => pointer on orbital in common list of orbitals
!---------------------------------------------------------------------
      Implicit none

      Integer :: ne     = 0     !  number of electrons
      Integer :: parity = 0     !  parity of states (+1,-1)

! ... description of 1 conf.w.function:

      Integer, parameter :: msh = 31 ! max. number of shells behind core

      Integer :: no, Jtotal, iconf, iterm
      Integer, dimension(msh) :: nn,kn,ln,jn,iq,in,       &
                                 Jshell,Vshell,Jintra, &
                                 np_symc,np_symt, np_orb, np,mp

      Integer :: no1, Jtotal1, iconf1, iterm1
      Integer, dimension(msh) :: nn1,kn1,ln1,jn1,iq1,in1,       &
                                 Jshell1,Vshell1,Jintra1, &
                                 np_symc1,np_symt1, np_orb1, np1,mp1
      Integer :: no2, Jtotal2, iconf2, iterm2
      Integer, dimension(msh) :: nn2,kn2,ln2,jn2,iq2,in2,       &
                                 Jshell2,Vshell2,Jintra2, &
                                 np_symc2,np_symt2, np_orb2, np2,mp2

! ... storing configuration as a character string:

      Integer, parameter :: mas = 9*msh+3
      Character(mas) :: CONFIG, SHELLJ, INTRAJ, AS,BS,CS
      Integer :: ia

! ... core:

      Integer,parameter :: mcore = 50
      Character(250) :: core, closed
      Integer :: ncore = 0
      Integer :: nn_core(mcore),k_core(mcore),l_core(mcore),j_core(mcore)
      Character(5) :: e_core(mcore)

! ... configuration list:

      Integer :: ncfg  = 0       !  current number of configurations
      Integer :: mcfg  = 0       !  max. number of configurations
      Integer :: icfg  = 50000   !  initial prediction of mcfg
      Integer :: jcfg  = 10      !  average number of shells
      Integer :: kcfg  = 0       !  max. dimension (mcfg*jcfg)
      Integer :: lcfg  = 0       !  last element

      Integer,    allocatable :: IP_state(:)
      Integer,    allocatable :: IS_term(:)
      Integer(2), allocatable :: IP_orb(:)

      Integer, allocatable :: IT_state1(:),IT_state2(:),IS_order(:)

      Integer :: J_min=-1, J_max=-1

! ... expansion coeficients:

      Real(8), allocatable :: WC(:)

! ... label assignments:

      Integer, parameter :: mlab = 64
      Character(mlab), allocatable :: LABEL(:)

! ... J-blocks:

      Integer :: njbl  = 0
      Integer, allocatable :: JJc (:),JTc1 (:),JTc2 (:),Jncfg (:), JTp(:)
      Integer :: njbl1 = 0
      Integer, allocatable :: JJ1c(:),JT1c1(:),JT1c2(:),Jncfg1(:), JTp1(:)
      Integer :: njbl2 = 0
      Integer, allocatable :: JJ2c(:),JT2c1(:),JT2c2(:),Jncfg2(:), JTp2(:)

      Real(8) :: memory_conf_jj = 0.d0

      End Module conf_jj


!======================================================================
      Subroutine alloc_cfg(m)
!======================================================================
!     allocate, deallocate or reallocate the arrays in module conf_jj
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer :: m,i
      Integer, allocatable :: iarr(:)
      Real(8), allocatable :: rarr(:)

      if(m.le.0) then
       if(allocated(ip_state)) Deallocate(ip_state,ip_orb,IS_term,WC)
       mcfg=0; ncfg = 0; lcfg = 0; ne = 0
      elseif(.not.allocated(ip_state)) then
       mcfg = m; kcfg = mcfg*jcfg
       Allocate(ip_state(mcfg),IS_term(mcfg),ip_orb(kcfg),WC(mcfg))
       ncfg = 0; lcfg = 0
      elseif(m.le.mcfg) then
       Return
      elseif(ncfg.le.0) then
       Deallocate (ip_state,ip_orb,IS_term,WC); ncfg = 0; lcfg = 0
       Allocate(ip_state(mcfg),IS_term(mcfg),ip_orb(kcfg),WC(mcfg))
      else
       mcfg=m; i=lcfg/ncfg+1; if(i.gt.jcfg) jcfg=i; kcfg=mcfg*jcfg
       Allocate(iarr(lcfg))
       iarr(1:ncfg)=IS_term(1:ncfg); Deallocate(IS_term)
       Allocate(IS_term(m)); IS_term=0; IS_term(1:ncfg)=iarr(1:ncfg)
       iarr(1:ncfg)=ip_state(1:ncfg); Deallocate(ip_state)
       Allocate(ip_state(m)); ip_state=0; ip_state(1:ncfg)=iarr(1:ncfg)
       iarr(1:lcfg)=ip_orb(1:lcfg); Deallocate(ip_orb)
       Allocate(ip_orb(kcfg)); ip_orb=0; ip_orb(1:lcfg)=iarr(1:lcfg)
       Deallocate(iarr)

       Allocate(rarr(ncfg))
       rarr(1:ncfg)=WC(1:ncfg); Deallocate(WC)
       Allocate(WC(m)); WC=0.d0; WC(1:ncfg)=rarr(1:ncfg)
       Deallocate(rarr)

      end if

      End Subroutine alloc_cfg


!======================================================================
      Integer Function Iadd_cfg_jj(job)
!======================================================================
!     add new or detect existing CAS in conf_jj list;
!     returns position of given CAS in the list
!     job  =  'add'     -  just add
!          =  'detect'  -  return -ic if exist
!          =   others   -  add if not exist
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Character(*), intent(in) :: job
      Integer :: i,ic,ip
      Integer, external :: Ifind_jjorb, Iadd_symc, Iadd_symt

      Iadd_cfg_jj = 0
      if(no.le.0) Return

      if(mcfg.eq.0) Call Alloc_cfg(icfg)
      Jtotal = Jintra(no)

!      Do i=2,no
!       if(nn(i-1).eq.nn(i).and.ln(i-1).eq.ln(i)) iq(i) = -iq(i)         ! ???
!      End do

      iconf = Iadd_symc(Jtotal,no,iq,kn)

!      iq = iabs(iq)                                                     ! ???

      iterm = Iadd_symt(iconf,no,Jshell,Vshell,Jintra)

      Do i = 1,no; np(i)=Ifind_jjorb(nn(i),kn(i),in(i),2); End do

! ... check if we already have such state:

      if(job.ne.'add') then
      Do ic = 1,ncfg
       if(IS_term(ic).ne.iterm) Cycle
       ip = ip_state(ic); Iadd_cfg_jj = ic; if(job.eq.'detect') Iadd_cfg_jj=-ic
       Do i = 1,no; ip=ip+1
        if(np(i).ne.IP_orb(ip)) then; Iadd_cfg_jj=0; Exit; end if
       End do
       if(Iadd_cfg_jj.ne.0) Return
      End do
      end if

      ncfg=ncfg+1
      if(ncfg.eq.mcfg.or.lcfg+no.gt.kcfg) Call Alloc_cfg(mcfg+icfg)
      IS_term(ncfg)=iterm
      ip_state(ncfg)=lcfg
      Do i=1,no; lcfg=lcfg+1; IP_orb(lcfg)=np(i); End do
      Iadd_cfg_jj = ncfg

      End Function Iadd_cfg_jj


!======================================================================
      Subroutine Get_cfg_jj(ic)
!======================================================================
!     extract the CAS "ic" from conf_jj list
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer :: ic, i,j,ip, ii

      iterm=IS_term(ic)
      Call Get_symt(iterm,iconf,no,Jshell,Vshell,Jintra)
      Call Get_symc(iconf,JTOTAL,no,nn,kn,ln,jn,iq,in)
      ip = ip_state(ic)
      Do i=1,no; ip=ip+1
       j = IP_orb(ip)
       Call Get_orb_jj(j,nn(i),kn(i),in(i),ii)
      End do

      End Subroutine Get_cfg_jj


!=======================================================================
      Integer Function Iort_conf_jj(kk)
!=======================================================================
!     define "strong" orthogonality between config.1 and config.2
!     based on the difference in kappa-values
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer :: i,j,k,kk

      Iort_conf_jj = 0
      np1=iq1; np2=iq2
      Do i=1,no1
       Do j=1,no2
        if(np2(j).eq.0) Cycle
        if(kn1(i).ne.kn2(j)) Cycle
        k=min(np1(i),np2(j)); np1(i)=np1(i)-k; np2(j)=np2(j)-k
        if(np1(i).eq.0) Exit
       End do
      End do
      k = SUM(np2(1:no2));  if(k.gt.kk) Iort_conf_jj = 1

      End Function Iort_conf_jj


!=======================================================================
      Subroutine Def_Jblocks
!=======================================================================
!     define the number of J-blocks in the conf_jj list, njbl,
!     and their attributes:
!     JJc   -  total angular momentum
!     Jncgf -  number of ACF in the block
!     JTc1  -  index of the first ACS in the block
!     JTc2  -  index of the lastt ACS in the block
!     JTp   -  parity of the block
!
!     remark: don't be confused - JJc and Jncfg first used as auxiliary
!             arrays to save the memory (re-write ?)
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer :: i,j

      Allocate(JJc(ncfg),Jncfg(ncfg))
      Call Get_cfg_jj(1)
      njbl=1; JJc(1) = Jtotal; Jncfg(1)=1

      Do i=2,ncfg
       Call Get_cfg_jj(i)
       if(Jtotal.eq.JJc(njbl)) then
        Jncfg(njbl)=i
       else
        njbl=njbl+1
        JJc(njbl)=Jtotal
        Jncfg(njbl)=i
       end if
      End do

      Allocate(JTc1(njbl),JTc2(njbl))
      JTc1(1:njbl)=JJc(1:njbl)
      JTc2(1:njbl)=Jncfg(1:njbl)
      Deallocate(JJc,Jncfg)
      Allocate(JJc(njbl),Jncfg(njbl))
      JJc=JTc1; JTc1(1)=1; Jncfg(1)=JTc2(1)
      Do i = 2,njbl
       JTc1(i) = JTc2(i-1)+1
       Jncfg(i) = JTc2(i) - JTc2(i-1)
      End do

      Allocate(JTp(njbl))
      Do i = 1,njbl
       Call Get_cfg_jj(JTc1(i))
       j=SUM(ln(1:no)*iq(1:no)); JTp(i)=(-1)**j
      End do

      End Subroutine Def_Jblocks


!======================================================================
      Integer Function no_ic_jj (ic)
!======================================================================
!     number of shells in state "ic"
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer, intent(in) :: ic
      Integer, external :: Get_iconf, Get_no

      iterm=IS_term(ic)
      iconf=Get_iconf(iterm)
      no_ic_jj = Get_no(iconf)

      End Function no_ic_jj


!======================================================================
      Integer Function jot_ic (ic)
!======================================================================
!     total momentum of state "ic"
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer, intent(in) :: ic
      Integer, external :: Get_iconf, Get_jot

      iterm=IS_term(ic)
      iconf=Get_iconf(iterm)
      jot_ic = Get_jot(iconf)

      End Function jot_ic

!======================================================================
      Subroutine Save_cfg_jj(i)
!======================================================================
!     save(restore) curent state in position i = (1|2)
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer, Intent(in) :: i

      Select Case(i)
      Case(1)
       no1=no; nn1=nn; kn1=kn; ln1=ln; jn1=jn; iq1=iq; in1=in
       Jshell1=Jshell; Vshell1=Vshell; Jintra1=Jintra
       Jtotal1=Jtotal; iconf1=iconf; iterm1=iterm
      Case(2)
       no2=no; nn2=nn; kn2=kn; ln2=ln; jn2=jn; iq2=iq; in2=in
       Jshell2=Jshell; Vshell2=Vshell; Jintra2=Jintra
       Jtotal2=Jtotal; iconf2=iconf; iterm2=iterm
      Case(-1)
       no=no1; nn=nn1; kn=kn1; ln=ln1; jn=jn1; iq=iq1; in=in1
       Jshell=Jshell1; Vshell=Vshell1; Jintra=Jintra1
       Jtotal=Jtotal1; iconf=iconf1; iterm=iterm1
      Case(-2)
       no=no2; nn=nn2; kn=kn2; ln=ln2; jn=jn2; iq=iq2; in=in2
       Jshell=Jshell2; Vshell=Vshell2; Jintra=Jintra2
       Jtotal=Jtotal2; iconf=iconf2; iterm=iterm2
      Case Default
       Stop 'Save_cfg_jj: ???'
      End Select

      End Subroutine Save_cfg_jj
