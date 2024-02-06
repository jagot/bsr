!===========================================================================
      Module channels
!===========================================================================
!     contains scattering channels information
!---------------------------------------------------------------------------
      Implicit none

      Integer :: nlsp =  0                !  number of partial waves

      Integer, allocatable :: ispar(:)    !  (2*S+1) for partial waves
      Integer, allocatable :: lpar(:)     !  total L
      Integer, allocatable :: ipar(:)     !  parity (+1 | -1)
      Integer, allocatable :: jpar(:)     !  (2*J+1)
      Integer, allocatable :: nch(:)      !  number of one-electron channels

      Integer, allocatable :: iptar(:,:)  !  pointer on the target state
      Integer, allocatable :: lch(:,:)    !  small l for given channel
      Integer, allocatable :: ipch(:,:)   !  pointer on the place in the common list of orbitals

      Integer, allocatable :: ipconf(:,:) !  pointer on the last configuration
      Integer, allocatable :: jkch(:,:)   !  k-number in jK-coupling - (2*J+1) ???

      CHARACTER(3), allocatable :: Tpar(:)  ! spectroscopic notation for partial wave
      CHARACTER(4), allocatable :: ELC(:,:) ! spectroscopic symbol for channel orbital

      Integer :: mch  =  0 !  maximum number of one-electron channels

! ... perturbers:

      Character(20), allocatable :: AFP(:) ! file-name for perturber
      Character(20), allocatable :: BFP(:) ! file-name for perturber

      Integer, allocatable :: ncp(:)  !  number of configurations in perturber
      Integer, allocatable :: nwp(:)  !  number of orbitals in perturber

      Integer, allocatable :: npert(:)   ! number of pertubers in given partial wave
      Integer, allocatable :: ipert(:)   ! base-pointer for given partial wave
      Integer, allocatable :: ippert(:)  ! pointer to the last configuration
      Integer :: mpert =0                ! max.dimension for ippert

      End Module channels


!=======================================================================
      Subroutine Allocate_channels
!=======================================================================

      Use channels

      if(nlsp.le.0) Stop ' Allocate_channels: nlsp <= 0 '
      if(mch.le.0) Stop ' Allocate_channels: mch <= 0 '

      Allocate(ispar(nlsp),lpar(nlsp),ipar(nlsp),jpar(nlsp),Tpar(nlsp), &
             nch(nlsp),iptar(nlsp,mch),lch(nlsp,mch),ipch(nlsp,mch),  &
             ELC(nlsp,mch), ipconf(nlsp,mch), jkch(nlsp,mch),         &
             AFP(nlsp),BFP(nlsp),ncp(nlsp),nwp(nlsp))
      lch = 0

      End subroutine Allocate_channels


!======================================================================
      Subroutine R_channels(nut)
!======================================================================
!     read channels information from file 'nut'
!----------------------------------------------------------------------
      Use channels

      Implicit none

      Integer, Intent(in) :: nut
      Character(20) :: AF
      Character(80) :: line
      Integer :: i,j,m,nu, Icheck_file, Ifind_position

      Call Read_ipar(nut,'max_ch',mch)
      if(mch.le.0) Stop 'R_channel: mch <= 0 '

      Call Read_ipar(nut,'nlsp',nlsp); read(nut,*)
      if(nlsp.le.0) Stop 'R_channel: nlsp <= 0 '

      Call Allocate_channels

      Do i = 1,nlsp
       read(nut,'(a)') line
       read(line,*) AF,lpar(i),ispar(i),ipar(i)
       AFP(i)=' '; ncp(i)=0; nwp(i)=0
       j = INDEX(line,'pert')
       if(j.gt.0) read(line(j:),*) AF,ncp(i),nwp(i)
       j=LEN_TRIM(AF); if(AF(j-1:j).eq.'.c') AF(j-1:j)='  '
       AFP(i)=AF
       jpar(i) = 0; if(ispar(i).eq.0) jpar(i) = lpar(i) + 1
      End do

      i = Ifind_position(nut,'channels:'); read(nut,*)

      Do i = 1,nlsp
       read(nut,*) line
       read(nut,*) AF,AF,AF,AF,nch(i),AF,AF,m,ncp(i)
       Do j = 1,nch(i)
        read(nut,*) ELC(i,j),lch(i,j),iptar(i,j),m,ipconf(i,j),jkch(i,j)
       End do
      End do

! ... pertubers:

      mpert = SUM(ncp)
      Allocate(npert(nlsp),ipert(nlsp),ippert(0:mpert))
      ipert = 0; npert = 0; ippert(0:mpert)=0
      if(mpert.eq.0) Return

      m = 0
      Do i = 1,nlsp
       ipert(i) = m
       if(ncp(i).eq.0) Cycle
       Call Find_free_unit(nu)
       AF = AFP(i); j=LEN_TRIM(AF); AF=AF(1:j)//'.c'
       if(Icheck_file(AF).eq.0) Cycle
       open(nu,file=AF)
       Call Read_ipar(nu,'npert',npert(i))
       if(npert(i).gt.0) then
        read(nu,*) ippert(m+1:m+npert(i)); m=m+npert(i)
       else
        npert(i)=ncp(i)
        Do j=1,npert(i); m=m+1; ippert(m)=j; End do
       end if
       Close(nu)
      End do

      End Subroutine R_channels


!======================================================================
      Subroutine Write_channels_LS (nut,met,max_nc,max_wf)
!======================================================================
!     write channel information in unit "nut"
!----------------------------------------------------------------------
      Use channels

      Integer, intent(in) :: nut,met,max_nc,max_wf

      write(nut,'(a,i4,5x,a)') &
           'nlsp  = ',nlsp,' !   number of partial waves'
      write(nut,'(80(''-''))')
      Do i = 1,nlsp
       write(nut,'(i3,3i4,6x,a40,2x,a10,2i5)') &
         i,lpar(i),ispar(i),ipar(i),AFP(i),BFP(i),ncp(i),nwp(i)
      End do
      write(nut,'(80(''-''))')

      if(met.eq.0) Return

      write(nut,'(a,i5)') 'channels:'
      write(nut,'(80(''-''))')

      Do ilsp = 1,nlsp
       write(nut,'(i3,a,i5.3,a,i6,a,2i10)') ilsp,'.',ilsp,' nch = ',nch(ilsp), &
        ' nc = ',ipconf(ilsp,nch(ilsp)),ncp(ilsp)
       Do i = 1,nch(ilsp)
        write(nut,'(a4,2x,3i6,2i8)') &
         ELC(ilsp,i),lch(ilsp,i),iptar(ilsp,i),i,ipconf(ilsp,i),jkch(ilsp,i)
       End do
       write(nut,'(80(''-''))')
      End do

      write(nut,'(a,i9)') 'max_ch =',mch
      write(nut,'(a,i9)') 'max_nc =',max_nc
      write(nut,'(a,i9)') 'max_wf =',max_wf
      write(nut,'(80(''-''))')

      End Subroutine Write_channels_LS
