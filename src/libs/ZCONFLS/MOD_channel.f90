!===========================================================================
      Module channel
!===========================================================================
!     define scattering channels for one partial wave
!---------------------------------------------------------------------------
      Implicit none

      Integer :: lpar =  0  ! total L for given patial wave
      Integer :: ispar=  0  ! (2*S+1) for given patial wave
      Integer :: ipar =  0  ! parity
      Integer :: jpar =  0  ! (2*J+1)
      Integer :: nch  =  0  ! number of one-electron channels
      Integer :: mch  =  0  ! max. number of channels
      Integer :: imch =256  ! initial prediction of mch
      Integer :: jmch =256  ! increment for mch

      Integer :: nLS  =  0  ! current dimension for LS arrays
      Integer :: mLS  =  0  ! max. dimension for LS-arrays

      Integer, allocatable :: iptar(:)  ! target pointer
      Integer, allocatable :: lch(:)    ! small l for given channel
      Integer, allocatable :: ipch(:)   ! place in the common list of orbitals
      Integer, allocatable :: ipconf(:) ! pointer to the last configuration
      Integer, allocatable :: jkch(:)   ! k-number in case of jK-coupling (2K+1)

      Integer, allocatable :: ip_LS(:)  ! pointer to the last configuration
      Integer, allocatable :: IL_LS(:)  ! pointer to the last configuration
      Integer, allocatable :: IS_LS(:)  ! pointer to the last configuration
      Integer, allocatable :: JL_LS(:)  ! pointer to the last configuration
      Integer, allocatable :: JS_LS(:)  ! pointer to the last configuration
      Real(8), allocatable :: CT_LS(:)  ! pointer to the last configuration

      Character(4), allocatable :: ELC(:) ! spectroscopic symbol

      Integer :: ncp	  ! number of configurations in perturber
      Integer :: nwp        ! number of orbitals in perturber
      Character(20) :: AFP  ! file-name for perturber

      Integer :: npert=0    ! number of configurations in perturber
      Integer :: mpert=0
      Integer :: ipert=640
      Integer, allocatable :: ippert(:)  ! pointer to the last configuration

      Character(3) :: Tpar

      End Module channel


!=======================================================================
      Subroutine Allocate_channel(m)
!=======================================================================
!     allocate arrays in module "channel"
!-----------------------------------------------------------------------
      Use channel

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: ia(:)
      Character(4), allocatable :: aa(:)

      if(m.le.0) then
       if(Allocated(iptar)) Deallocate(iptar,lch,ipch,ELC,ipconf,jkch,ip_LS)
       mch = 0; nch=0
      elseif(.not.allocated(iptar)) then
       mch = m
       Allocate(iptar(mch),ipconf(mch),lch(mch),ipch(mch), &
                ELC(mch),jkch(mch),ip_LS(0:mch))
                ip_LS = 0
      elseif(m.le.mch) then
       Return
      elseif(nch.gt.0) then
       mch = m
       Allocate(ia(nch))
       ia=iptar(1:nch); Deallocate(iptar)
       Allocate(iptar(mch)); iptar(1:nch)=ia
       ia=ipconf(1:nch); Deallocate(ipconf)
       Allocate(ipconf(mch)); ipconf(1:nch)=ia
       ia=lch(1:nch); Deallocate(lch)
       Allocate(lch(mch)); lch(1:nch)=ia
       ia=ipch(1:nch); Deallocate(ipch)
       Allocate(ipch(mch)); ipch(1:nch)=ia
       ia=jkch(1:nch); Deallocate(jkch)
       Allocate(jkch(mch)); jkch(1:nch)=ia
       ia=ip_LS(1:nch); Deallocate(ip_LS)
       Allocate(ip_LS(0:mch)); ip_LS(1:nch)=ia; ip_LS(0)=0
       Deallocate(ia)
       Allocate(aa(nch))
       aa=ELC(1:nch); Deallocate(ELC)
       Allocate(ELC(mch)); ELC(1:nch)=aa
       Deallocate(aa)
      end if

      End Subroutine Allocate_channel


!=======================================================================
      Subroutine Allocate_LS_arrays(m)
!=======================================================================
      Use channel

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: ia(:)
      Real(8), allocatable :: ra(:)

      if(m.le.0) then
       if(Allocated(IL_LS)) Deallocate(IL_LS,IS_LS,JL_LS,JS_LS,CT_LS)
       mLS = 0
      elseif(.not.allocated(IL_LS)) then
       mLS = m
       Allocate(IL_LS(mLS),IS_LS(mLS),JL_LS(mLS),JS_LS(mLS),CT_LS(mLS))
      elseif(m.le.mLS) then
       Return
      elseif(nLS.gt.0) then
       mLS = m
       Allocate(ia(mLS))
       ia=IL_LS(1:nLS); Deallocate(IL_LS)
       Allocate(IL_LS(mLS)); IL_LS(1:nLS)=ia
       ia=IS_LS(1:nLS); Deallocate(IS_LS)
       Allocate(IS_LS(mLS)); IS_LS(1:nLS)=ia
       ia=JL_LS(1:nLS); Deallocate(JL_LS)
       Allocate(JL_LS(mLS)); JL_LS(1:nLS)=ia
       ia=JS_LS(1:nLS); Deallocate(JS_LS)
       Allocate(JS_LS(mLS)); JS_LS(1:nLS)=ia
       Deallocate(ia)
       Allocate(ra(nLS))
       ra=CT_LS(1:nLS); Deallocate(CT_LS)
       Allocate(CT_LS(mLS)); CT_LS(1:nLS)=ra
       Deallocate(ra)
      end if

      End Subroutine Allocate_LS_arrays


!=======================================================================
      Subroutine Allocate_pert(m)
!=======================================================================
!     allocate "pertuber" arrays in module "channel"
!-----------------------------------------------------------------------
      Use channel

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: ia(:)

      if(m.le.0) then
       if(Allocated(ippert)) Deallocate(ippert)
       mpert = 0; npert=0
      elseif(.not.allocated(ippert)) then
       mpert = m
       Allocate(ippert(0:mpert)); ippert=0
      elseif(m.le.mpert) then
       Return
      else
       mpert = m
       Allocate(ia(0:npert))
       ia=ippert(0:npert); Deallocate(ippert)
       Allocate(ippert(0:mpert)); ippert(0:npert)=ia(0:npert)
       Deallocate(ia)
      end if

      End Subroutine Allocate_pert


!======================================================================
      Subroutine R_channel(nut,klsp)
!======================================================================
!     read from file 'nut' information for channel klsp
!----------------------------------------------------------------------
      Use channel

      Implicit none
      Integer, intent(in) :: nut,klsp
      Character(20) :: AF
      Character(80) :: line
      Integer :: i,j, nlsp, ncfg, nu, Ifind_position

      Call Read_ipar(nut,'nlsp',nlsp); read(nut,*)
      if(nlsp.le.0) Stop 'R_channel: nlsp <= 0 '
      if(klsp.gt.nlsp) Stop 'R_channel: klsp > nlsp '

      Do i = 1,klsp
       read(nut,'(a)') line
       read(line,*) AF,lpar,ispar,ipar
       AFP=' '; ncp=0; nwp=0
       j = INDEX(line,'pert')
       if(j.gt.0) read(line(j:),*) AFP,ncp,nwp
       jpar = -1; if(ispar.eq.0) jpar = lpar + 1
      End do

      i = Ifind_position(nut,'channels:')
      read(nut,*)
      read(nut,*)

      Call Allocate_channel(0)
      Do
       read(nut,'(a)') line;
       if(line(4:4).ne.'.') Cycle
       read(line(1:3),*) i;  if(i.ne.klsp) Cycle
       read(line,*) AF,AF,AF,AF,nch,AF,AF,ncfg,ncp
       i = nch; Call Allocate_channel(i)
       Do i = 1,nch
        read(nut,*) ELC(i),lch(i),iptar(i),j,ipconf(i),jkch(i)
       End do
       Exit
      End do

      Call Allocate_pert(0)
      if(ncp.le.0) Return

      Call Find_free_unit(nu)
      i = LEN_TRIM(AFP); AF=AFP(1:i)//'.c'
      open(nu,file=AF)
      Call R_pert(nu)
      Close(nu)

      End Subroutine R_channel


!======================================================================
      Subroutine R_pert(nu)
!======================================================================
!     read from unit 'nu' information about perturbers
!----------------------------------------------------------------------
      Use channel

      Implicit none
      Integer, Intent(in) :: nu
      Integer :: i

      Call Allocate_pert(0)
      Call Read_ipar(nu,'npert',npert)
      if(npert.gt.0) then
       Call Allocate_pert(ipert+npert)
       read(nu,*) ippert(1:npert)
      elseif(ncp.gt.0) then
       npert=ncp
       Call Allocate_pert(ipert+npert)
       Do i=0,npert; ippert(i)=i; End do
      end if

      End Subroutine R_pert


!======================================================================
      Integer Function Ifind_channel(ic)
!======================================================================
!     find channel (or perturber) for given configuration 'ic'
!----------------------------------------------------------------------
      Use channel

      Implicit none
      Integer, Intent(in) :: ic
      Integer :: ich

      if(ic.le.0) Call Stop_mpi (0,0,'Ifind_channel: ic < 0')
      Ifind_channel = 1
      Do ich = 1,nch
       if(ic.gt.ipconf(ich)) Ifind_channel=ich+1
      End do

      if(Ifind_channel.le.nch) Return

      Do ich = 1,npert
       if(ic.gt.ippert(ich)) Ifind_channel=nch+ich+1
      End do

      if(Ifind_channel.gt.nch+npert) then
       write(*,*) 'npert=',npert
       write(*,*)  ippert(0:npert)
       Call Stop_mpi (0,0,'Ifind_channel: ic is to big')
      end if

      End Function Ifind_channel


!======================================================================
      Subroutine Find_channel_ic(ich,ic1,ic2)
!======================================================================
!     find configurations for given channel (or perturber)
!----------------------------------------------------------------------
      Use channel

      Implicit none
      Integer, intent(in)  :: ich
      Integer, intent(out) :: ic1,ic2
      Integer :: i

      if(ich.le.nch) then
       i = ich
       ic1=1; if(i.gt.1) ic1=ipconf(i-1)+1
       ic2=ipconf(i)
      else
       i = ich-nch
       ic1=ippert(i-1)+1
       ic2=ippert(i)
      end if

      End Subroutine Find_channel_ic

