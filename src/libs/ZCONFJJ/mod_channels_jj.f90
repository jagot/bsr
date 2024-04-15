!===========================================================================
      Module channels_jj
!===========================================================================
!     define scattering channels in jj-coupling
!---------------------------------------------------------------------------
      Implicit none

      Integer :: nlsp = 0    !   number of partial waves
      Integer :: mch  = 0    !   max. number of one-electron channels

      Integer, allocatable :: jpar(:) !  2*J - value
      Integer, allocatable :: ipar(:) !  parity (+1 | -1)
      Integer, allocatable :: nch(:)  !  number of one-electron channels

      Integer, allocatable :: iptar(:,:)
      Integer, allocatable :: kch(:,:)
      Integer, allocatable :: lch(:,:)
      Integer, allocatable :: jch(:,:)
      Integer, allocatable :: ipch(:,:)
      Integer, allocatable :: ipconf(:,:)

      !  iptar  -  pointer on the target state
      !  kch    -  kappa l for given channel
      !  lch    -  small l for given channel
      !  jch    -  2j-value for given channel
      !  ipch   -  pointer on the place in the common list of orbitals
      !  ipconf -  pointer on the last configuration

      Character(5), allocatable :: ELC(:,:)

      !  ELC     -  spectroscopic symbol for given channel orbital

      Character(40), allocatable :: AFP(:)
      Character(10), allocatable :: BFP(:)
      Integer, allocatable :: ncp(:)
      Integer, allocatable :: nwp(:)

      !  AFP     -  file-name for perturber (given)
      !  BFP     -  file-name for perturber (internal)
      !  ncp     -  number of configurations in perturber
      !  nwp     -  number of orbitals in perturber

      Integer, allocatable :: npert(:)   ! number of pertubers in given partial wave
      Integer, allocatable :: ipert(:)   ! base-pointer for given partial wave
      Integer, allocatable :: ippert(:)  ! pointer to the last configuration
      Integer :: mpert =0                ! max.dimension for ippert

      Integer :: max_nc = 0    !   max. number of configurations
      Integer :: max_wf = 0    !   max. number of one-electron w.f.

      End Module channels_jj


!=======================================================================
      Subroutine Alloc_channels_jj
!=======================================================================
!     just allocate arrays in module channels_jj
!---------------------------------------------------------------------
      Use channels_jj

      if(nlsp.le.0) Stop ' Alloc_channels_jj: nlsp <= 0 '
      if(mch .le.0) Stop ' Alloc_channels_jj: mch  <= 0 '

      Allocate(jpar(nlsp),ipar(nlsp),nch(nlsp), &
             iptar(nlsp,mch),kch(nlsp,mch),lch(nlsp,mch),jch(nlsp,mch), &
             ipch(nlsp,mch), ELC(nlsp,mch), ipconf(nlsp,mch), &
             AFP(nlsp),BFP(nlsp),ncp(nlsp),nwp(nlsp))

      End subroutine Alloc_channels_jj


!======================================================================
      Subroutine Read_channels_jj(nut)
!======================================================================
!     read channels information from unit 'nut'
!----------------------------------------------------------------------
      USE channels_jj

      Integer, intent(in) :: nut
      Character(20) :: AF
      Character(80) :: line
      Integer :: i,j,k
      Integer, external :: Ifind_position

      Call Read_ipar(nut,'max_ch',mch)
      if(mch.le.0) Stop 'R_channel: mch <= 0 '

      Call Read_ipar(nut,'nlsp',nlsp); read(nut,*)
      if(nlsp.le.0) Stop 'R_channel: nlsp <= 0 '

      Call Alloc_channels_jj

      AFP = ' '; BFP = ' '; ncp=0; nwp=0
      Do i = 1,nlsp
       read(nut,'(a)') line
       read(line,*) AF,jpar(i),ipar(i)
       if(LEN_TRIM(line).gt.12) read(line(13:),*) AFP(i),BFP(i),ncp(i),nwp(i)
      End do

      i = Ifind_position(nut,'channels:'); read(nut,*)

      Do i = 1,nlsp
       read(nut,*)
       read(nut,'(a)') line
       j = INDEX(line,'nch =') + 5
       read(line(j:),*) nch(i)
       read(nut,*)
       Do j = 1,nch(i)
        read(nut,'(a)') line
        ELC(i,j)=line(7:11)
        read(line(12:),*) kch(i,j),iptar(i,j),ipconf(i,j)
        lch(i,j) = l_kappa(kch(i,j))
        jch(i,j) = j_kappa(kch(i,j))
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
       AF = BFP(i); j=LEN_TRIM(AF); AF=AF(1:j)//'.c'
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

      End Subroutine Read_channels_jj


!======================================================================
      Subroutine Write_channels_jj (nut,met)
!======================================================================
!     write channel information in unit "nut"
!----------------------------------------------------------------------
      Use channels_jj

      Integer, intent(in) :: nut,met

      write(nut,'(a,i4,5x,a)') &
           'nlsp  = ',nlsp,' !   number of partial waves'
      write(nut,'(80(''-''))')

      Do i = 1,nlsp
       if(ncp(i).eq.0) then
        write(nut,'(i3,a,2i4,6x,a40,2x,a10,2i5)') &
         i,'.',jpar(i),ipar(i)
       else
        write(nut,'(i3,a,2i4,6x,a40,2x,a10,2i5)') &
         i,'.',jpar(i),ipar(i),AFP(i),BFP(i),ncp(i),nwp(i)
       end if
      End do
      write(nut,'(80(''-''))')

      if(met.eq.0) Return

      write(nut,'(a,i5)') 'channels: '
      write(nut,'(80(''-''))')

      Do ilsp = 1,nlsp
       write(nut,'(i3,a,a,i6)')  ilsp,'.',' nch =',nch(ilsp)
       write(nut,'(80(''-''))')
       Do i = 1,nch(ilsp)
        write(nut,'(i3,a,2x,a,2x,2i6,i8)') &
         i,'.',ELC(ilsp,i),kch(ilsp,i),iptar(ilsp,i),ipconf(ilsp,i)
       End do
       write(nut,'(80(''-''))')
      End do

      write(nut,'(a,i9)') 'max_ch =',mch
      write(nut,'(a,i9)') 'max_nc =',max_nc
      write(nut,'(a,i9)') 'max_wf =',max_wf
      write(nut,'(80(''-''))')

      End Subroutine Write_channels_jj


!======================================================================
      Integer Function Ifind_pert_jj(ilsp,ic)
!======================================================================
!     find perturber for given configuration 'ic'
!----------------------------------------------------------------------
      Use channels_jj

      Implicit none
      Integer, intent(in) :: ilsp,ic
      Integer :: i, ii, jc

      Ifind_pert_jj=ic;    if(ilsp.le.0) Return

      Ifind_pert_jj=0;     if(ic.le.ipconf(ilsp,nch(ilsp))) Return

      jc=ic-ipconf(ilsp,nch(ilsp)); ii = ipert(ilsp)

      Ifind_pert_jj = -1
      Do i = 1,npert(ilsp)
       if(jc.gt.ippert(ii+i)) Cycle
       Ifind_pert_jj = i; Exit
      End do

      if(Ifind_pert_jj.eq.-1) Stop 'Ifind_pert_jj: ic is to big'

      End Function Ifind_pert_jj
