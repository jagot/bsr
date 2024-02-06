!======================================================================
      Module target_ion
!======================================================================
!     contains description of target states;
!     it is just a copy of module "target" in case when we have two
!     scattering systems to be described
!----------------------------------------------------------------------
      Implicit none

      Integer :: ntarg	 !  number of target states
      Integer :: nphys	 !  number of physical states
      Integer :: nelc    !  number of atomic electrons
      Integer :: nz      !  atomic number
      Integer :: nct	 !  total number of target config.s
      Integer :: nwt     !  total number of taget orbitals

      Integer, allocatable :: istarg(:) !  (2*S+1) for target i
      Integer, allocatable :: ltarg(:)  !  total L
      Integer, allocatable :: iptarg(:) !  parity (+-1)
      Integer, allocatable :: jtarg(:)  !  (2*J+1)
      REAL(8), allocatable :: etarg(:)  !  target energy in au

      Integer, allocatable :: nctarg(:) !  number of target configurations
      Integer, allocatable :: nwtarg(:) !  number of new orbitals
      Integer, allocatable :: ictarg(:) !  pointer to target i in conf.list

      CHARACTER(20), allocatable :: AFT(:) ! given file-names for target states
      CHARACTER(20), allocatable :: BFT(:) ! original file-names for target states

      CHARACTER(2) :: COUPLING = 'LS' ! LS, JK or JJ coupling mode

      End Module target_ion


!======================================================================
      Subroutine allocate_target_ion(m)
!======================================================================
!     allocate (deallocate) space in module "target"
!----------------------------------------------------------------------
      Use target_ion

      Implicit none
      Integer, intent(in) :: m

      if(m.le.0) then
       if(allocated(istarg)) Deallocate(istarg,ltarg,iptarg,jtarg, &
                nctarg, nwtarg, ictarg, etarg, AFT, BFT)
       ntarg = 0
      else
       if(allocated(istarg)) Deallocate(istarg,ltarg,iptarg,jtarg, &
                nctarg, nwtarg, ictarg, etarg, AFT, BFT)
       ntarg = m
       ALLOCATE(istarg(ntarg),  ltarg(ntarg), iptarg(ntarg), &
                nctarg(ntarg), nwtarg(ntarg), ictarg(ntarg), &
                etarg(ntarg), AFT(ntarg),BFT(ntarg), jtarg(ntarg) )
      end if

      nct = 0
      nwt = 0

      End Subroutine allocate_target_ion


!======================================================================
      Subroutine R_target_ion(nut)
!======================================================================
!     read from file 'nut' target information
!----------------------------------------------------------------------
      Use target_ion

      Implicit none
      Integer, intent(in) :: nut
      Character(20) :: AF
      Integer :: i

      Call Read_ipar(nut,'nelc',nelc)
      Call Read_ipar(nut,'nz',nz)
      Call Read_ipar(nut,'nphys',nphys)
      Call Read_apar(nut,'coupling',coupling)
      Call Read_ipar(nut,'ntarg',ntarg)

      if(ntarg.le.0) Stop 'R_targ: ntarg <= 0 '
      i = ntarg; Call Allocate_target_ion(i)

      Call Read_ipar(nut,'nwt',nwt)
      Call Read_ipar(nut,'ntarg',ntarg)

      nct = 0
      read(nut,*)
      Do i=1,ntarg
       read(nut,*) BFT(i),AFT(i),ltarg(i),istarg(i),iptarg(i),etarg(i), &
                   nctarg(i),nwtarg(i)
       nct=nct+nctarg(i)
!       nwt=nwt+nwtarg(i)
       ictarg(i) = nct
       jtarg(i)=0; if(istarg(i).eq.0) jtarg(i)=ltarg(i)+1
      End do

      End Subroutine R_target_ion

