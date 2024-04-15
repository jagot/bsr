!======================================================================
      Module channel_jj
!======================================================================
! ... defines scattering channels for one partial wave
!----------------------------------------------------------------------
      Implicit none

      Integer :: jpar =  0  !  2*J
      Integer :: ipar =  0  !  parity
      Integer :: nch  =  0  !  number of one-electron channels
      Integer :: mch  =  0  !  max. number of channels
      Integer :: imch =512  !  initial prediction(or incriment) for mch

      Integer, allocatable :: kch(:)      ! kappa value for given channel
      Integer, allocatable :: lch(:)      ! l-value for given channel
      Integer, allocatable :: jjch(:)     ! l-value for given channel
      Integer, allocatable :: iptar(:)    ! associated target state
      Integer, allocatable :: ipconf(:)   ! last configuration for given channel
      Integer, allocatable :: ipch(:)     ! additional pointer
      Character(5), allocatable :: ELC(:) ! channel spectroscopic symbol

      Integer :: ncp	      !  number of configurations in perturber
      Integer :: nwp            !  number of orbitals in perturber
      Character(80) :: AFP,BFP  ! file-name for perturber

      Integer :: npert=0        ! number of configurations in perturber
      Integer :: mpert=0
      Integer :: ipert=64
      Integer, allocatable :: ippert(:)  ! pointer to the last configuration

      Character(4) :: Tpar !  some given notation for given partial wave

      End Module channel_jj


!======================================================================
      Subroutine Alloc_channel_jj(m)
!======================================================================
! ... allocate, deallocate or reallocate arrays in module channel_jj
!---------------------------------------------------------------------
      Use channel_jj

      Implicit none
      Integer :: m,i
      Integer, allocatable :: iarr(:)
      Character(5), allocatable :: arr(:)

      if(m.le.0) then
        if(allocated(iptar)) Deallocate(iptar,kch,lch,jjch,ipconf,ipch,ELC)
        mch = 0
      elseif(.not.allocated(iptar)) then
        mch = m
        Allocate(iptar(m),kch(m),lch(m),jjch(m),ipconf(m),ipch(m),ELC(m))
      elseif(m.le.mch) then
        Return
      elseif(nch.eq.0) then
        Deallocate(iptar,kch,lch,jjch,ipconf,ipch,ELC)
        mch = m
        Allocate(iptar(m),kch(m),lch(m),jjch(m),ipconf(m),ipch(m),ELC(m))
      else
        Allocate(iarr(m))
        iarr(1:nch)=iptar(1:nch); Deallocate(iptar)
        Allocate(iptar(m)); iptar(1:nch)=iarr(1:nch)
        iarr(1:nch)=kch(1:nch); Deallocate(kch)
        Allocate(kch(m)); kch(1:nch)=iarr(1:nch)
        iarr(1:nch)=lch(1:nch); Deallocate(lch)
        Allocate(lch(m)); lch(1:nch)=iarr(1:nch)
        iarr(1:nch)=jjch(1:nch); Deallocate(jjch)
        Allocate(jjch(m)); jjch(1:nch)=iarr(1:nch)
        iarr(1:nch)=ipconf(1:nch); Deallocate(ipconf)
        Allocate(ipconf(m)); ipconf(1:nch)=iarr(1:nch)
        iarr(1:nch)=ipch(1:nch); Deallocate(ipch)
        Allocate(ipch(m)); ipch(1:nch)=iarr(1:nch)
        Deallocate(iarr)
        Allocate(arr(m))
        arr(1:nch)=ELC(1:nch); Deallocate(ELC)
        Allocate(ELC(m)); ELC(1:nch)=arr(1:nch)
        Deallocate(arr)
        mch = m
      end if

      End Subroutine Alloc_channel_jj


!======================================================================
      Subroutine Read_channel_jj(nut,klsp)
!======================================================================
! ... read from file 'nut' information for channel klsp
!----------------------------------------------------------------------
      USE channel_jj

      Implicit none
      Integer, Intent(in) :: nut,klsp
      Character(20) :: AF
      Character(80) :: line
      Integer :: i,j, nlsp, ncfg, nu
      Integer, external :: l_kappa, j_kappa, Ifind_position

      Call Read_ipar(nut,'nlsp',nlsp); read(nut,*)
      if(nlsp.le.0) Stop 'R_channel: nlsp <= 0 '
      if(klsp.gt.nlsp) Stop 'R_channel: klsp > nlsp '

      Do i = 1,klsp; read(nut,'(a)') line;  End do
      read(line,*) Tpar,jpar,ipar
      AFP = ' '; BFP = ' '; ncp=0; nwp=0
      if(LEN_TRIM(line).gt.12) read(line(13:),*) AFP,BFP,ncp,nwp

      i = Ifind_position(nut,'channels:'); read(nut,*)
      Do
       read(nut,*); read(nut,'(a)') line; read(nut,*)
       read(line(1:3),*) i
       j = INDEX(line,'nch =') + 5
       read(line(j:),*) nch
       if(i.ne.klsp) then
        Do i=1,nch; read(nut,'(a)') line; End do; Cycle
       end if
       Call Alloc_channel_jj(0); Call Alloc_channel_jj(nch)
       Do i = 1,nch; read(nut,'(a)') line
        ELC(i)=line(7:11); read(line(12:),*) kch(i),iptar(i),ipconf(i)
        lch(i) = l_kappa(kch(i))
        jjch(i) = j_kappa(kch(i))
       End do
       Exit
      End do

      Call Allocate_pert(0)
      if(ncp.le.0) Return
      Call Find_free_unit(nu)
      i = LEN_TRIM(BFP); AF=BFP(1:i)//'.c'
      open(nu,file=AF)
      Call R_pert(nu)

      Close(nu)

      End Subroutine Read_channel_jj


!======================================================================
      Subroutine R_pert(nu)
!======================================================================
!     read from unit 'nu' information about perturbers
!----------------------------------------------------------------------
      Use channel_jj

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
!       ippert = ippert + ipconf(nch)     ???

      End Subroutine R_pert

!=======================================================================
      Subroutine Allocate_pert(m)
!=======================================================================
!     allocate "pertuber" arrays in module "channel"
!-----------------------------------------------------------------------
      Use channel_jj

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
      Integer Function Ifind_channel_jj(ic)
!======================================================================
!     find channel (or perturber) for given configuration 'ic'
!----------------------------------------------------------------------
      Use channel_jj

      Implicit none
      Integer, intent(in) :: ic
      Integer :: ich

      if(ic.le.0) Stop 'Ifind_channel_jj: ic <= 0'
      Ifind_channel_jj = 1
      Do ich = 1,nch
       if(ic.gt.ipconf(ich)) Ifind_channel_jj=ich+1
      End do

      if(Ifind_channel_jj.le.nch) Return

      Do ich = 1,npert
       if(ic.gt.ippert(ich)) Ifind_channel_jj=nch+ich+1
      End do

      if(Ifind_channel_jj.gt.nch+npert) then
       write(*,*) 'npert=',npert
       write(*,*)  ippert(0:npert)
       Stop 'Ifind_channel_jj: configuration index, ic, is to big'
      end if

      End Function Ifind_channel_jj
