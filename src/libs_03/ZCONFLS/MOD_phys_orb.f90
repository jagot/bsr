!=================================================================
      Module phys_orb_LS
!=================================================================
!     contains the list of physical orbitals for all target states
!-----------------------------------------------------------------
      Implicit none

      Integer :: nphys_orb = 0           !  number of phys.orbitals

      Integer, allocatable :: ip_tar(:)  !  target pointer
      Integer, allocatable :: ip_phy(:)  !  orbital pointer
      Integer, allocatable :: ip_sub(:)  !  substitution pointer
      Real(8), allocatable :: S_orb (:)  !  population

      Integer :: nphys_sub = 0           !  number of sub.orbitals
      Integer, allocatable :: jp_sub(:)  !  index of sub. orbitals
      Integer, allocatable :: jp_tar(:)  !  index of sub. orbitals

      Integer :: npert_sub = 0           !  number of sub.orbitals
      Integer, allocatable :: np_phy(:)  !  orbital pointer
      Integer, allocatable :: np_sub(:)  !  substitution pointer
      Real(8), allocatable :: SP_orb(:)  !  population

      End Module phys_orb_LS


!================================================================
      Subroutine read_sub_orb_LS(nu,ntarg)
!================================================================
!     read orbitals information from unit 'nu'
!----------------------------------------------------------------
      Use phys_orb_LS

      Implicit none
      Integer, intent(in) :: nu, ntarg
      Character(80) :: A
      Integer :: i,j,n,l,k,kap,jot, it,ip,i1,i2
      Integer, external :: Ifind_nlk, Ipointer, Ifind_position

! ... define nphys_orb

      nphys_orb = 0
      if(Ifind_position(nu,'target').eq.0) Return

      Do i=1,ntarg
       read(nu,*)
       Do
        read(nu,'(a)') A
        if(A(1:1).eq.'*') Exit
        nphys_orb=nphys_orb+1
       End do
      End do

      if(allocated(ip_tar)) Deallocate(ip_tar,ip_phy,ip_sub,S_orb)
      Allocate(ip_tar(ntarg),ip_phy(nphys_orb),ip_sub(nphys_orb), &
               S_orb(nphys_orb))

! ... define physical and corr. substitution orbitals orbitals:

      i=Ifind_position(nu,'target'); j=0
      Do i=1,ntarg
       read(nu,*)
       Do
        read(nu,'(a)') A
        if(A(1:1).eq.'*') Exit
        j=j+1
        Call EL4_nlk(A(1:4),n,l,k)
        ip_phy(j) = Ifind_nlk(n,l,k,2)
        Call EL4_nlk(A(18:21),n,l,k)
        ip_sub(j) = Ifind_nlk(n,l,k,2)
        read(A(11:),*) S_orb(j)
       End do
       ip_tar(i) = j
      End do

! ... define list of all different substitution orbitals:

      if(allocated(jp_sub)) Deallocate(jp_sub)
      Allocate(jp_sub(nphys_orb))
      nphys_sub=1; jp_sub(1)=ip_sub(1)
      Do i=2,nphys_orb
       if(Ipointer(nphys_sub,jp_sub,ip_sub(i)).ne.0) Cycle
       nphys_sub=nphys_sub+1
       jp_sub(nphys_sub)=ip_sub(i)
      End do

! ... define target states where the substitution orbital
!     first appears:

      if(allocated(jp_tar)) Deallocate(jp_tar)
      Allocate(jp_tar(nphys_orb))
      jp_tar = 0
      Do it=1,ntarg
       i1=1; if(it.gt.1) i1=ip_tar(it-1)+1; i2=ip_tar(it)
       Do i=i1,i2
        ip = Ipointer(nphys_sub,jp_sub,ip_sub(i))
        if(jp_tar(ip).eq.0) jp_tar(ip)=it
       End do
      End do

      End Subroutine read_sub_orb_LS


!================================================================
      Subroutine read_sub_pert_LS(nu,ipert)
!================================================================
      Use phys_orb_LS

      Implicit none
      Integer, Intent(in) :: nu, ipert
      Character(80) :: A
      Integer :: i,j,n,l,k,kap,jot,jpert
      Integer, Allocatable :: iarr(:)
      Integer, External :: Ifind_nlk, Ipointer

! ... define npert_sub

      npert_sub = 0
      rewind(nu)
      Do
       read(nu,'(a)',end=10) A
       if(A(1:6).ne.'pertub') Cycle
       read(A(7:),*) jpert
       if(jpert.ne.ipert) Cycle
       Do
        read(nu,'(a)') A
        if(A(1:1).eq.'*') Exit
        npert_sub = npert_sub + 1
       End Do
       Exit
      End Do
      if(npert_sub.eq.0) Return
      if(Allocated(np_phy)) Deallocate(np_phy,np_sub,SP_orb)
      Allocate(np_phy(npert_sub),np_sub(npert_sub),SP_orb(npert_sub))

! ... define physical and corr. substitution orbitals orbitals:

      rewind(nu)
      j = 0
      Do
       read(nu,'(a)',end=10) A
       if(A(1:6).ne.'pertub') Cycle
       read(A(7:),*) jpert
       if(jpert.ne.ipert) Cycle
       Do
        read(nu,'(a)') A
        if(A(1:1).eq.'*') Exit
        j=j+1
        Call EL4_nlk(A(1:4),n,l,k)
        np_phy(j) = Ifind_nlk(n,l,k,2)
        Call EL4_nlk(A(18:21),n,l,k)
        np_sub(j) = Ifind_nlk(n,l,k,2)
        read(A(11:),*) SP_orb(j)
       End do
       Exit
      End Do
   10 Continue

      End Subroutine read_sub_pert_LS

