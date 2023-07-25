!======================================================================
      MODULE phys_orb_jj
!======================================================================
!     contains the list of physical and substitution orbitals 
!     for each target state: "ip_phy" and "ip_sub" 
!     jp_sub - just total list of substitution orbitals 
!----------------------------------------------------------------------
      Implicit none
   
      Integer :: nphys_orb = 0           !  number of phys.orbitals

      Integer, allocatable :: ip_tar(:)  !  target pointer 
      Integer, allocatable :: ip_phy(:)  !  orbital pointer 
      Integer, allocatable :: ip_sub(:)  !  substitution pointer 
      Real(8), allocatable :: S_orb (:)  !  occupation

      Integer :: nphys_sub = 0           !  number of all sub.orbitals
      Integer, allocatable :: jp_sub(:)  !  index of sub. orbitals 

      Integer :: npert_sub = 0           !  number of sub.orbitals for pertuber
      Integer, allocatable :: np_phy(:)  !  orbital pointer 
      Integer, allocatable :: np_sub(:)  !  substitution pointer 
      Real(8), allocatable :: SP_orb (:) !  occupation

      END MODULE phys_orb_jj


!======================================================================
      Subroutine read_sub_orb_jj(nu,ntarg)
!======================================================================
!     read substitution orbitals from unit "nu" for "ntarg" states
!----------------------------------------------------------------------
      USE phys_orb_jj

      Implicit none
      Integer, intent(in) :: nu, ntarg
      Character(80) :: A
      Integer :: i,j,n,l,k,kap,jot
      Integer, external :: Ifind_jjorb, Ipointer, Ifind_position

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
        Call EL_nljk(A(1:5),n,kap,l,jot,k)
        ip_phy(j) = Ifind_jjorb(n,kap,k,2)
        Call EL_nljk(A(19:23),n,kap,l,jot,k)
        ip_sub(j) = Ifind_jjorb(n,kap,k,2)
        read(A(11:),*) S_orb(j)
       End do
       ip_tar(i) = j
      End do

! ... define list of all different substitution orbitals:

      if(allocated(jp_sub)) deallocate(jp_sub)
      allocate(jp_sub(nphys_orb))
      nphys_sub=1; jp_sub(1)=ip_sub(1)
      Do i=2,nphys_orb
       if(Ipointer(nphys_sub,jp_sub,ip_sub(i)).ne.0) Cycle
       nphys_sub=nphys_sub+1
       jp_sub(nphys_sub)=ip_sub(i)
      End do

      End Subroutine read_sub_orb_jj


!======================================================================
      Subroutine read_sub_pert_jj(nu,ipert)
!======================================================================
!     read substitution orbitals from unit "nu" for perturber "ipert"
!----------------------------------------------------------------------
      Use phys_orb_jj

      Implicit none
      Integer, intent(in) :: nu, ipert
      Character(80) :: A
      Integer :: i,j,n,l,k,kap,jot,jpert
      Integer, external :: Ifind_jjorb, Ipointer

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
        Call EL_nljk(A(1:5),n,kap,l,jot,k)
        np_phy(j) = Ifind_jjorb(n,kap,k,2)
        Call EL_nljk(A(19:23),n,kap,l,jot,k)
        np_sub(j) = Ifind_jjorb(n,kap,k,2)
        read(A(11:),*) SP_orb(j)
       End do
       Exit
      End Do

   10 Continue

      End Subroutine read_sub_pert_jj

