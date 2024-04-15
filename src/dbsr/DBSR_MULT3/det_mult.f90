!======================================================================
      Subroutine Det_mult
!======================================================================
!     creates the common list of orbital symmetries for two input
!     determinants and call the subroutines for calculations of
!     m.e. between possible combinations of nj-orbitals
!----------------------------------------------------------------------
      Use dbsr_mult
      Use nljm_orbitals;  Use conf_jj, only: ne

      Implicit none
      Integer :: i,j, i1,j1, k, k1,k2
      Integer, external :: Isort
!----------------------------------------------------------------------
!                              creat common list of orbital symmetries:
      ksym1=1; ksym2=1

      Nsym=0; k1=0;  k2=0; kz1=0; kz2=0

! ... exhaust the 1-st configuration:

      Do i = 1,ne

       if(ksym1(i).eq.0) Cycle
       Nsym = Nsym + 1
       Lsym(Nsym)=Lsym1(i); Msym(Nsym)=Msym1(i); Jsym(Nsym)=Jsym1(i)
       k1=k1+1; IPsym1(Nsym)=k1; Isym1(k1)=i; ksym1(i)=0

! ... check for the same orbitals rest the 1-st configuration:

       Do j = i+1,ne
        if(ksym1(j).eq.0) Cycle
        if(Lsym(Nsym).ne.Lsym1(j)) Cycle
        if(Msym(Nsym).ne.Msym1(j)) Cycle
        if(Jsym(Nsym).ne.Jsym1(j)) Cycle
        k1=k1+1; IPsym1(Nsym)=k1; Isym1(k1)=j; ksym1(j)=0
       End do

!       Jdet=Isym1(1:ne);  kz1 = Isort(ne,Jdet) ???

! ... check for the same orbitals the 2-nd configuration:

       IPsym2(Nsym)=k2
       Do j = 1,ne
        if(ksym2(j).eq.0) Cycle
        if(Lsym(Nsym).ne.Lsym2(j)) Cycle
        if(Msym(Nsym).ne.Msym2(j)) Cycle
        if(Jsym(Nsym).ne.Jsym2(j)) Cycle
        k2=k2+1; IPsym2(Nsym)=k2; Isym2(k2)=j; ksym2(j)=0
       End do

      End do

      if(k1.ne.ne) Stop 'Det_me: k1 <> ne '

! ... exhaust the 2-st configuration:

      Do i = 1,ne
       if(ksym2(i).eq.0) Cycle
       Nsym = Nsym + 1
       Lsym(Nsym)=Lsym2(i); Msym(Nsym)=Msym2(i); Jsym(Nsym)=Jsym2(i)
       k2=k2+1; IPsym2(Nsym)=k2; Isym2(k2)=i; ksym2(i)=0
       IPsym1(Nsym)=k1

! ... check for the same orbitals rest of 2-st configuration:

       Do j = i+1,ne
        if(ksym2(j).eq.0) Cycle
        if(Lsym(Nsym).ne.Lsym2(j)) Cycle
        if(Msym(Nsym).ne.Msym2(j)) Cycle
        if(Jsym(Nsym).ne.Jsym2(j)) Cycle
        k2=k2+1; IPsym2(Nsym)=k2; Isym2(k2)=j; ksym2(j)=0
       End do

      End do

      if(k2.ne.ne) Stop 'Det_breit: k2 <> ne '

      Jdet=Isym1(1:ne);  kz1 = Isort(ne,Jdet)
      Jdet=Isym2(1:ne);  kz2 = Isort(ne,Jdet)

!----------------------------------------------------------------------
!                              define the number of different orbitals:
      Ksym1(1)=ipsym1(1)
      Ksym2(1)=ipsym2(1)
      Do i = 2,NSYM
       Ksym1(i)=ipsym1(i)-ipsym1(i-1)
       Ksym2(i)=ipsym2(i)-ipsym2(i-1)
      End do

! ... how much different symmetries:

      k = 0
      Do i = 1,NSYM
       N1(i) = KSYM1(i)-KSYM2(i)
       N2(i) = KSYM2(i)-KSYM1(i)
       if(N1(i).gt.0) k = k + N1(i)
      End do

      if(k.gt.1) Return

      Select case (k)
!---------------------------------------------------------------------
!                                                         k = 1  case:
      Case(1)

       if(kpol.eq.0) Return

       Do i=1,NSYM; if(N1(i).le.0) Cycle; i1=i; Exit;  End do
       Do i=1,NSYM; if(N2(i).le.0) Cycle; j1=i; Exit;  End do

       Call ZNO_001(i1,j1)

!---------------------------------------------------------------------
!                                                         k = 0  case:
      Case(0)

      if(kpol.eq.0) then
       Call ZNO_overlap
      else
       Call ZNO_000
      end if

      End Select

      End Subroutine DET_mult


!====================================================================
      Subroutine ZNO_001(is,js)
!====================================================================
! ... angular part of electric multipole operator between two det.w.f
!--------------------------------------------------------------------
      Use dbsr_mult; Use nljm_orbitals

      Implicit none
      Real(8) :: CA,CB
      Integer, intent(in) :: is,js
      Integer :: i,j,i1,i2,j1,j2,k1,k2,idf,int,kz
      Integer, external :: Idet_fact, Incode_mult

      Call me_jj(Lsym(is),Jsym(is),Msym(is), &
                 Lsym(js),Jsym(js),Msym(js),CA,CB)

      if(CA.eq.0.d0) Return

      i1 = 1; if(is.gt.1) i1=IPsym1(is-1)+1; i2=IPsym1(is)
      j1 = 1; if(js.gt.1) j1=IPsym2(js-1)+1; j2=IPsym2(js)
      Do i=i1,i2; k1=nnsym1(Isym1(i))
      Do j=j1,j2; k2=nnsym2(Isym2(j))

       idf = Idet_fact(i,0,j,0);  kz = (-1)**(kz1+kz2+i+j)

       Select case(ktype)
        Case('E')
         int = Incode_mult(1,k1,k2);  Call Iadd_zoef(CA*kz,int,idf)
        Case('M')
         int = Incode_mult(2,k1,k2);  Call Iadd_zoef(CA*kz,int,idf)
       End Select

      End do; End do

      End Subroutine ZNO_001


!====================================================================
      Subroutine ZNO_000
!====================================================================
!    angular part of electric multipole operator between two det.w.f
!--------------------------------------------------------------------
      Use dbsr_mult; Use nljm_orbitals

      Implicit none
      Real(8) :: CA,CB
      Integer :: i,j,i1,i2,k,k1,k2,is,idf,int,kz
      Integer, external :: Idet_fact, Incode_mult

      Do is = 1,NSYM

      Call me_jj(Lsym(is),Jsym(is),Msym(is), &
                 Lsym(is),Jsym(is),Msym(is),CA,CB)

       if(CA.eq.0.d0) Cycle

       i1 = 1; if(is.gt.1) i1=IPsym1(is-1)+1; i2=IPsym1(is)
       Do i=i1,i2; k=Isym1(i); k1=nnsym1(k)
       Do j=i1,i2; k=Isym2(j); k2=nnsym2(k)

       idf = Idet_fact(i,0,j,0); kz = (-1)**(kz1+kz2+i+j)

       Select case(ktype)
        Case('E')
         int = Incode_mult(1,k1,k2);  Call Iadd_zoef(CA*kz,int,idf)
        Case('M')
         int = Incode_mult(2,k1,k2);  Call Iadd_zoef(CA*kz,int,idf)
       End Select

       End do; End do

      End do

      End Subroutine ZNO_000


!====================================================================
      Subroutine ZNO_overlap
!====================================================================
!     computes overlap integral between two determinants
!
!     Calls: Idet_fact, Iadd_zoef, Incode_mult
!--------------------------------------------------------------------
      Use nljm_orbitals

      Implicit none
      Integer :: idf,int
      Real(8) :: C
      Integer, external :: Idet_fact, Incode_mult

      C = (-1)**(kz1+kz2)
      idf = Idet_fact (0,0,0,0)
      int = Incode_mult (0,1,1)

      Call Iadd_zoef (C,int,idf)

      End Subroutine ZNO_overlap


