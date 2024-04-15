!======================================================================
      Module target_jj
!======================================================================
!     contains description of target states
!----------------------------------------------------------------------
      Implicit none

      Character(80) :: TITLE

      Integer :: ntarg       ! number of target states
      Integer :: nelc	     ! number of atomic electrons
      Integer :: nz          ! atomic number
      Integer :: nct	     ! total number of target config.s
      Integer :: nwt         ! total number of taget orbitals

      Real(8), allocatable :: etarg(:)   ! target energy in au
      Integer, allocatable :: jtarg(:)   ! 2J-value
      Integer, allocatable :: ptarg(:)   ! parity (+1|-1)
      Integer, allocatable :: nctarg(:)  ! number of target configurations
      Integer, allocatable :: nwtarg(:)  ! number of new orbitals
      Integer, allocatable :: ictarg(:)  ! pointer to target in conf. list

      Character(30), allocatable :: AFT(:)  ! given file-names for target states
      Character(10), allocatable :: BFT(:)  ! original file-names for target states

      End Module target_jj


!======================================================================
      Subroutine alloc_target_jj(m)
!======================================================================
!     allocate (deallocate) space in module "target_jj"
!----------------------------------------------------------------------
      Use target_jj

      Implicit none
      Integer, Intent(in) :: m

      if(m.le.0) then
       if(allocated(jtarg)) Deallocate(jtarg,ptarg, &
                            nctarg, nwtarg, ictarg, etarg, AFT,BFT)
       ntarg = 0
      else
       if(allocated(jtarg)) Deallocate(jtarg,ptarg, &
                            nctarg, nwtarg, ictarg, etarg, AFT,BFT)
       ntarg = m
       ALLOCATE(jtarg(ntarg),ptarg(ntarg),etarg(ntarg), &
                nctarg(ntarg), nwtarg(ntarg), ictarg(ntarg), &
                AFT(ntarg),BFT(ntarg))
      end if
      nct = 0
      nwt = 0

      End Subroutine alloc_target_jj


!======================================================================
      Subroutine Read_target_jj(nut)
!======================================================================
!     read target information from unit 'nut'
!----------------------------------------------------------------------
      Use target_jj

      Implicit none
      Integer, intent(in) :: nut
      Integer :: i

      rewind(nut)
      read(nut,'(a)') title

      Call Read_ipar(nut,'nelc',nelc)
      Call Read_ipar(nut,'nz',nz)
      Call Read_ipar(nut,'ntarg',ntarg)

      if(ntarg.le.0) Stop 'R_targ: ntarg <= 0 '
      i=ntarg; Call Alloc_target_jj(i)

      Call Read_ipar(nut,'ntarg',ntarg)
      read(nut,*)
      nct = 0
      Do i=1,ntarg
       read(nut,*) AFT(i),BFT(i),jtarg(i),ptarg(i),etarg(i), &
                   nctarg(i),nwtarg(i)
       nct=nct+nctarg(i)
       ictarg(i) = nct
      End do
      nwt = 0; Call Read_ipar(nut,'nwt',nwt)

      End Subroutine Read_target_jj


!======================================================================
      Subroutine Write_target_jj(nut)
!======================================================================
!     write target information to file 'nut'
!----------------------------------------------------------------------
      Use target_jj

      Implicit none
      Integer, Intent(in) :: nut
      Integer :: i

      rewind(nut)
      write(nut,'(a)') TRIM(title)
      write(nut,'(80(''-''))')
      write(nut,'(a,i4,5x,a)') &
                'nz    = ',nz,   ' !   nuclear charge'
      write(nut,'(a,i4,5x,a)') &
                'nelc  = ',nelc, ' !   number of electrons'
      write(nut,'(80(''-''))')
      write(nut,'(a,i4,5x,a)') &
                'ntarg = ',ntarg,' !   number of target states'
      write(nut,'(80(''-''))')
      Do i=1,ntarg
       write(nut,'(a30,2x,a10,2x,2i4,f18.8,2i5)') AFT(i),BFT(i), &
        jtarg(i),ptarg(i),etarg(i),nctarg(i),nwtarg(i)
      End do
      write(nut,'(80(''-''))')
      write(nut,'(a,i7)') 'nct =',nct
      write(nut,'(a,i7)') 'nwt =',nwt
      write(nut,'(80(''-''))')

      End Subroutine Write_target_jj
