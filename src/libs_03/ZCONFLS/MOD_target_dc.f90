!======================================================================
      Module target_dc
!======================================================================
!     contains description of target for double continuum case
!----------------------------------------------------------------------
      Implicit none

      Integer :: ntarg	    !  number of target states
      Integer :: ntar1
      Integer :: ntar2
      Integer :: nelc	      !  number of atomic electrons
      Integer :: nz         !  atomic number
      Integer :: nct	       !  total number of target config.s
      Integer :: nwt        !  total number of taget orbitals

! ... istarg(i) - (2*S+1) for target i
! ... ltarg (i) - total L
! ... iptarg(i) - parity (+-1)

      Integer, allocatable :: istarg(:), ltarg(:), iptarg(:)

! ... etarg(i)  - target energy in au

      Real(8), allocatable :: etarg(:)

! ... nctarg(i)  - number of target configurations
! ... nwtarg(i)  - number of new orbitals
! ... ictarg(i)  - pointer to target i in conf.list

      Integer, allocatable :: nctarg(:), nwtarg(:), ictarg(:)

! ... AFT    - file-names for target states

      Character(20), allocatable :: AFT(:), BFT(:)

      End Module target_dc


!======================================================================
      Subroutine alloc_target_dc(m)
!======================================================================
!     allocate (deallocate) space in module "target"
!----------------------------------------------------------------------
      Use target_dc

      Implicit none
      Integer, intent(in) :: m

      if(m.le.0) then
       if(allocated(istarg)) Deallocate(istarg,ltarg,iptarg, &
                nctarg, nwtarg, ictarg, etarg, AFT, BFT)
       ntarg=0; ntar1=0; ntar2=0
      else
       if(allocated(istarg)) Deallocate(istarg,ltarg,iptarg, &
                nctarg, nwtarg, ictarg, etarg, AFT, BFT)
       ntarg = m
       ALLOCATE(istarg(ntarg),  ltarg(ntarg), iptarg(ntarg), &
                nctarg(ntarg), nwtarg(ntarg), ictarg(0:ntarg), &
                etarg(ntarg), AFT(ntarg),BFT(ntarg) )
      end if

      nct = 0
      nwt = 0

      End Subroutine alloc_target_dc


!======================================================================
      Subroutine Read_targ_dc(nut)
!======================================================================
!     read from file 'nut' target information
!----------------------------------------------------------------------
      Use target_dc

      Implicit none
      Integer, intent(in) :: nut
      Character(20) :: AF
      Integer :: i

      Call Read_ipar(nut,'nelc',nelc)
      Call Read_ipar(nut,'nz',nz)
      Call Read_ipar(nut,'ntar1',ntar1)
      Call Read_ipar(nut,'ntar2',ntar2)
      ntarg = ntar1 + ntar2
      if(ntarg.le.0) Stop 'R_targ_dc: ntarg <= 0 '
      Call Alloc_target_dc(ntar1+ntar2)

      nct=0; nwt=0; ictarg=0
      Call Read_ipar(nut,'ntar1',ntar1)
      read(nut,*)
      Do i=1,ntar1
       read(nut,*) BFT(i),AFT(i),ltarg(i),istarg(i),iptarg(i),etarg(i), &
                   nctarg(i),nwtarg(i)
       nct=nct+nctarg(i)
       nwt=nwt+nwtarg(i)
       ictarg(i) = nct
      End do

      Call Read_ipar(nut,'ntar2',ntar2)
      read(nut,*)
      Do i=ntar1+1,ntarg
       read(nut,*) BFT(i),AFT(i),ltarg(i),istarg(i),iptarg(i),etarg(i), &
                   nctarg(i),nwtarg(i)
       nct=nct+nctarg(i)
       nwt=nwt+nwtarg(i)
       ictarg(i) = nct
      End do

      End Subroutine Read_targ_dc

