!======================================================================
      Module dets
!======================================================================
!     contains information about overlap determinants {<i|j>}
!----------------------------------------------------------------------
!
! KPD(i)  - size of i-th determinant
!
! IPD(i)  - pointer (ip) of i-th determinants in the NPD array
!
! NPD(ip+1:ip+kd) - contains information about orbitals in overlap
!                   determinant as  ii * idet + jj, 
!                   where ii - pointer on row orbitals
!                         jj - pointer on column orbitals
!----------------------------------------------------------------------     
!
! KPF(i)  - number (kd) of determinants in i-th overlap factor
!
! IPF(i)  - pointer (ip) on the i-th overlap factor in the NPF array
!
! NPF(ip+1:ip+kd) - i-th overlap factors as a list of pointers 
!                   on individual determinants (ip) and its power (iext): 
!                   npf(i) = ip * idef  +  iext
!----------------------------------------------------------------------
! parameters idet,idef should agree with program bsr_breit
!----------------------------------------------------------------------      
      Implicit none

      Integer :: ndet  =   0   !  number of determinants
      Integer :: kdet  =   0   !  sum of all det. dimensions      
      Integer :: jmdet =   0
      Integer, parameter :: idet  = 2**15 !  pack basis 
	
      Integer, allocatable :: KPD(:), IPD(:), NPD(:), JPD(:)

      Integer :: ndef   =   0  ! number of overlap factors 
      Integer :: kdef   =   0  ! sum of all overlap factor dimensions  
      Integer :: jmdef  =   0
      Integer, parameter :: idef   =  16  ! pack basis 
      
      Integer, allocatable :: KPF(:), IPF(:), NPF(:), JPF(:)

      End Module dets


!======================================================================
      Subroutine Load_dets (nu)
!======================================================================
!     Containes the overlap determinants for all symmetries
!----------------------------------------------------------------------
      Use dets 

      Implicit none
      Integer :: nu

      read(nu) ndet,kdet,jmdet

      if(ndet.eq.0) Return

      if(allocated(kpd)) Deallocate (KPD,IPD,NPD,JPD)
      Allocate(KPD(ndet),IPD(ndet),NPD(kdet),JPD(ndet))

      read(nu) kpd
      read(nu) ipd
      read(nu) npd
      read(nu) jpd

      read(nu) ndef,kdef,jmdef
      if(ndef.eq.0) Return

      if(allocated(kpf)) Deallocate (KPF,IPF,NPF,JPF)
      Allocate(KPF(ndef),IPF(ndef),NPF(kdef),JPF(ndef))

      read(nu) kpf
      read(nu) ipf
      read(nu) npf
      read(nu) jpf

      End Subroutine Load_dets 

