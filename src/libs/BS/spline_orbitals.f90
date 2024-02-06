!=======================================================================
      Module spline_orbitals
!=======================================================================
!     contains description of atomic orbitals
!     and their B-spline representation
!-----------------------------------------------------------------------
      Implicit none

! ... LIST OF ONE-ELECTRON ORBITALS

      Integer :: mbf = 0               ! max. number of orbitals
      Integer :: nbf = 0               ! current number of orbitals
      Integer :: ibf = 512             ! initial prediction of mbf
      Integer :: jbf = 512             ! incriment for mbf

      Integer, allocatable :: nbs(:)   !  n-values
      Integer, allocatable :: lbs(:)   !  l-values
      Integer, allocatable :: kbs(:)   !  set numbers
      Integer, allocatable :: mbs(:)   !  number of splines
      Integer, allocatable :: iech(:)  !  additional pointer

      CHARACTER(4), allocatable :: ebs(:) ! spectroscopic notation

! ... B-spline expansion coefficients:

      Real(8), allocatable :: PBS(:,:)

! ... convolution with B-overlaps:

      Real(8), allocatable :: QBS(:,:)

! ... memory requirements (in 4b words):

      Integer :: m_borb = 0
      Real(8) :: memory_BS_orbitals = 0.d0

! ... overlaps between bound and continuum orbitals:
! ... (with <.|i> vectors if any):

      Integer :: nv_ch = 0
      Real(8), allocatable :: V_ch(:,:)
      Integer, allocatable :: i_ch(:),j_ch(:)
      Integer, allocatable :: my_channel(:)
      Integer, parameter   :: ib = 2**15
      Real(8) :: eps_v = 1.d-10

      End Module spline_orbitals


!=======================================================================
      Subroutine allocate_bsorb(m)
!=======================================================================
!     This program allocates (deallocates) space for list of atomic
!     orbitals or reallocates it if necessary
!-----------------------------------------------------------------------
      Use spline_orbitals
      Use spline_param

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: iarr1(:)
      Real(8), allocatable :: rarr2(:,:)
      Character(4), allocatable :: abs(:)

      if(m.le.0) then

       if(Allocated(NBS)) &
         Deallocate (NBS,LBS,KBS,MBS,iech,EBS,PBS,QBS)
         nbf = 0;  mbf = 0

      elseif(m.gt.mbf.and.nbf.eq.0) then

       if(Allocated(nbs)) &
          Deallocate (nbs,lbs,kbs,mbs,iech,ebs,PBS,QBS)

       mbf = m
       Allocate(nbs(mbf),lbs(mbf),kbs(mbf),ebs(mbf),mbs(1:mbf), &
                iech(1:mbf),PBS(1:ns,1:mbf), QBS(1:ns,1:mbf))
       nbs = 0; lbs = 0; kbs = 0; ebs = '****'; mbs = 0; iech = 0
       PBS = 0.d0; QBS = 0.d0

      elseif(nbf.gt.0.and.m.gt.mbf) then

       Allocate(iarr1(nbf))
       iarr1(1:nbf)=nbs(1:nbf); Deallocate(nbs);  Allocate(nbs(m))
       nbs=0;  nbs(1:nbf)=iarr1(1:nbf)
       iarr1(1:nbf)=lbs(1:nbf); Deallocate(lbs);  Allocate(lbs(m))
       lbs=0;  lbs(1:nbf)=iarr1(1:nbf)
       iarr1(1:nbf)=kbs(1:nbf); Deallocate(kbs);  Allocate(kbs(m))
       kbs=0;  kbs(1:nbf)=iarr1(1:nbf)
       iarr1(1:nbf)=mbs(1:nbf); Deallocate(mbs);  Allocate(mbs(m))
       mbs=0;  mbs(1:nbf)=iarr1(1:nbf)
       iarr1(1:nbf)=iech(1:nbf); Deallocate(iech);  Allocate(iech(m))
       iech=0; iech(1:nbf)=iarr1(1:nbf)
       Deallocate(iarr1)

       Allocate(abs(nbf))
       abs(1:nbf)=ebs(1:nbf); Deallocate(ebs);  Allocate(ebs(m))
       ebs='****'; ebs(1:nbf)=abs(1:nbf); Deallocate(abs)

       Allocate(rarr2(ns,nbf))
       rarr2(1:ns,1:nbf)=pbs(1:ns,1:nbf); Deallocate(pbs)
       Allocate(pbs(ns,m))
       pbs=0.d0; pbs(1:ns,1:nbf)=rarr2(1:ns,1:nbf)
       rarr2(1:ns,1:nbf)=qbs(1:ns,1:nbf); Deallocate(qbs)
       Allocate(qbs(ns,m))
       qbs=0.d0; qbs(1:ns,1:nbf)=rarr2(1:ns,1:nbf)
       Deallocate(rarr2)

       mbf=m;  ! write(*,*) 'realoc_BS_orb: mbf=',mbf

      end if

      m_borb = 6*mbf + mbf*mbf + 2 * (2*ns*mbf + mbf*mbf) + 6

      End Subroutine allocate_bsorb


!=======================================================================
      Integer function Ifind_bsorb(n,l,k)
!=======================================================================
! ... find the orbital nlk in the list, returns 0 if no such orbitals
!-----------------------------------------------------------------------
      USE spline_orbitals

      Implicit none
      Integer, intent(in) :: n,l,k
      Integer :: i

      Ifind_bsorb=0

      Do i=1,nbf
       if(n.eq.nbs(i).and.l.eq.lbs(i).and.k.eq.kbs(i)) then
        Ifind_bsorb = i
        Return
       end if
      End do

      END function Ifind_bsorb


!=======================================================================
      Integer function Jfind_bsorb(n,l,k)
!=======================================================================
! ... find the orbital nlk in the list, stops if no such orbitals
!-----------------------------------------------------------------------
      USE spline_orbitals

      Implicit none
      Integer, INTENT(in) :: n,l,k
      Integer, External :: Ifind_bsorb
      Integer :: i

      i =  Ifind_bsorb(n,l,k)
      if(i.eq.0) then
       Write(*,*) ' Jfind_bsorb: can not find the orbital: N,L,K=',n,l,k
       Stop
      else
       Jfind_bsorb = i
      end if

      END function Jfind_bsorb


!=======================================================================
      Integer function Iadd_bsorb(n,l,k)
!=======================================================================
!     adds the orbital to the list if new one
!-----------------------------------------------------------------------
      USE spline_orbitals

      Implicit none
      Integer, intent(in) :: n,l,k
      Integer, external :: Ifind_bsorb
      CHARACTER(4), EXTERNAL :: ELF4
      Integer :: i

      i =  Ifind_bsorb(n,l,k)

      if(i.eq.0) then

       if(nbf+1.gt.mbf) Call Allocate_bsorb(mbf+jbf)
       nbf = nbf + 1
       nbs(nbf) = n
       lbs(nbf) = l
       kbs(nbf) = k
       ebs(nbf) = ELF4(n,l,k)
       Iadd_bsorb = nbf

      else

       Iadd_bsorb = i

      end if

      End function Iadd_bsorb


