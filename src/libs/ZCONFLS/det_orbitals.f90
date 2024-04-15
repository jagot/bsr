!======================================================================
      Subroutine DET_orbitals1
!======================================================================
!     creates the common list of orbital symmetries for two input
!     determinants and call the subroutines for calculations of
!     coefficients between possible combinations of tese symmetries
!
!     Isym  - new disribution of orbitals
!     IPsym - pointer on last same orbutal in the new list
!----------------------------------------------------------------------
      Use spin_orbitals
      Use conf_LS,         only: ne

      Implicit none
      Integer :: i, j, k1
      Integer, external :: Isort

      NSYM = 0
!----------------------------------------------------------------------
! ... exzaust the 1-st configuration:

      ksym1=1; k1=0

      Do i = 1,ne
       if(ksym1(i).eq.0) Cycle
       Nsym = Nsym + 1
       Lsym(Nsym)=Lsym1(i); Msym(Nsym)=Msym1(i); Ssym(Nsym)=Ssym1(i)
       k1=k1+1; IPsym1(Nsym)=k1; Isym1(k1)=i; ksym1(i)=0

! ... check for the same orbitals rest the 1-st configuration:

       Do j = i+1,ne
        if(ksym1(j).eq.0) Cycle
        if(Lsym(Nsym).ne.Lsym1(j)) Cycle
        if(Msym(Nsym).ne.Msym1(j)) Cycle
        if(Ssym(Nsym).ne.Ssym1(j)) Cycle
        k1=k1+1; IPsym1(Nsym)=k1; Isym1(k1)=j; ksym1(j)=0
       End do

      End do

      if(k1.ne.ne) Stop 'Det_breit: k1 <> ne '
      IN=Isym1(1:ne);  kz1 = Isort(ne,IN)
      NSYM1 = NSYM

      End Subroutine DET_orbitals1


!======================================================================
      Subroutine DET_orbitals2
!======================================================================
!     creates the common list of orbital symmetries for two input
!     determinants and call the subroutines for calculations of
!     coefficients between possible combinations of symmetries
!
!     Isym - new disribution (perutatiom) of orbitals
!     IPsym - pointer on last same orbutal in new list
!----------------------------------------------------------------------
      Use spin_orbitals
      Use conf_LS, only: ne, joper

      Implicit none
      Integer :: i,i1,i2, j,j1,j2, k,k2, N1(2*ne), N2(2*ne)
      Integer, external :: Isort

      NSYM = nsym1
!----------------------------------------------------------------------
! ... check for the same orbitals in the 2-nd configuration:

      ksym2 = 1;  k2 = 0; IPsym2=0

      Do i = 1,NSYM
       Do j = 1,ne
        if(ksym2(j).eq.0) Cycle
        if(Lsym(i).ne.Lsym2(j)) Cycle
        if(Msym(i).ne.Msym2(j)) Cycle
        if(Ssym(i).ne.Ssym2(j)) Cycle
        k2=k2+1; IPsym2(i)=k2; Isym2(k2)=j; ksym2(j)=0
       End do
      End do

! ... exzaust the 2-st configuration:

      Do i = 1,ne
       if(ksym2(i).eq.0) Cycle
       Nsym = Nsym + 1
       Lsym(Nsym)=Lsym2(i); Msym(Nsym)=Msym2(i); Ssym(Nsym)=Ssym2(i)
       k2=k2+1; IPsym2(Nsym)=k2; Isym2(k2)=i; ksym2(i)=0

! ... check for the same orbitals rest of 2-st configuration:

       Do j = i+1,ne
        if(ksym2(j).eq.0) Cycle
        if(Lsym(Nsym).ne.Lsym2(j)) Cycle
        if(Msym(Nsym).ne.Msym2(j)) Cycle
        if(Ssym(Nsym).ne.Ssym2(j)) Cycle
        k2=k2+1; IPsym2(Nsym)=k2; Isym2(k2)=j; ksym2(j)=0
       End do
      End do

      if(k2.ne.ne) Stop 'Det_breit: k2 <> ne '
      IN=Isym2(1:ne);  kz2 = Isort(ne,IN)


      if(Nsym.gt.Nsym1) IPsym1(nsym1+1:nsym)=IPsym1(nsym1)
      Do i=2,nsym
       if(IPsym2(i).eq.0) IPsym2(i)=IPsym2(i-1)
      End do

!----------------------------------------------------------------------
!                              define the number of different orbitals:
      Ksym1(1)=ipsym1(1)
      Ksym2(1)=ipsym2(1)
      Do i = 2,NSYM
       Ksym1(i)=ipsym1(i)-ipsym1(i-1)
       Ksym2(i)=ipsym2(i)-ipsym2(i-1)
      End do

! ... how much symmetries differ:

      k = 0
      Do i = 1,NSYM
       N1(i) = KSYM1(i)-KSYM2(i)
       N2(i) = KSYM2(i)-KSYM1(i)
       if(N1(i).gt.0) k = k + N1(i)
      End do

      if(k.gt.2) Return

!---------------------------------------------------------------------
      Select case(k)

      Case(2)

       if(joper(3)+joper(5)+joper(6)+joper(7).eq.0) Return

       Do i=1,NSYM
        if(N1(i).le.0) Cycle; i1=i; N1(i)=N1(i)-1; Exit
       End do
       Do i=i1,NSYM
        if(N1(i).le.0) Cycle; i2=i; Exit
       End do
       Do i=1,NSYM
        if(N2(i).le.0) Cycle; j1=i; N2(i)=N2(i)-1; Exit
       End do
       Do i=j1,NSYM
        if(N2(i).le.0) Cycle; j2=i; Exit
       End do

       Call Zno_2ee(i1,i2,j1,j2)

      Case(1)

       if(joper(3)+joper(5)+joper(6)+joper(7).eq.0) Return

       Do i=1,NSYM; if(N1(i).le.0) Cycle; i1 = i; Exit; End do
       Do i=1,NSYM; if(N2(i).le.0) Cycle; j1 = i; Exit; End do

       Do i = 1,Nsym
        if(Ksym1(i).eq.0) Cycle
        if(i.eq.i1.and.Ksym1(i).le.1) Cycle
        if(i.eq.j1.and.Ksym2(i).le.1) Cycle

        if(i.le.i1.and.i.le.j1)  then
          Call Zno_2ee(i,i1,i,j1)
        elseif(i.gt.i1.and.i.le.j1) then
          Call Zno_2ee(i1,i,i,j1)
        elseif(i.gt.i1.and.i.gt.j1) then
          Call Zno_2ee(i1,i,j1,i)
        elseif(i.le.i1.and.i.gt.j1) then
          Call Zno_2ee(i,i1,j1,i)
        end if

       End do

      Case(0)

       if(joper(1).gt.0)           Call ZNO_0ee
       if(joper(2)+joper(4).gt.0)  Call ZNO_1ee

       if(joper(3)+joper(5)+joper(6)+joper(7).eq.0) Return

       Do i = 1,Nsym
        Do j = i,Nsym
         if(i.eq.j.and.Ksym1(i).le.1) Cycle
         Call Zno_2ee(i,j,i,j)
        End do
       End do

      End Select

      End Subroutine DET_orbitals2


!====================================================================
      Subroutine ZNO_0ee
!====================================================================
!     computes overlap integral between two determinants
!     Calls: Idet_fact, Incode_int, Iadd_zoef
!--------------------------------------------------------------------
      Use spin_orbitals,  only: kz1,kz2

      Implicit none
      Integer :: idf,int
      Real(8) :: C
      Integer, external :: Idet_fact, Incode_int

      C = (-1)**(kz1+kz2)
      idf = Idet_fact (0,0,0,0)
      int = Incode_int (11,0,1,1,1,1)

      Call Iadd_zoef (C,int,idf)

      End Subroutine ZNO_0ee

!====================================================================
      Subroutine ZNO_1ee
!====================================================================
!     angular part of one-electron operator between two det.w.f
!     Calls: Idet_fact, Incode_int, Iadd_zoef.
!--------------------------------------------------------------------
      Use conf_LS,       only: joper
      Use spin_orbitals, only: nsym, kz1,kz2,Msym,Ssym,IPsym1, &
                               Isym1,Isym2,nnsym1,nnsym2
      Implicit none
      Integer :: i,j,i1,i2,k,k1,k2,is,idf,int
      Real(8) :: C,CFF
      Integer, external :: Idet_fact, Incode_int

      Do is = 1,NSYM
       i1 = 1; if(is.gt.1) i1=IPsym1(is-1)+1; i2=IPsym1(is)

       Do i=i1,i2; k=Isym1(i); k1=nnsym1(k)
       Do j=i1,i2; k=Isym2(j); k2=nnsym2(k)

       idf = Idet_fact(i,0,j,0)

       C=(-1)**(kz1+kz2+i+j)

       if(joper(2).gt.0) then                          ! L-integrals
        int = Incode_int (6,0,k1,k1,k2,k2)
        CFF=-0.5;  Call Iadd_zoef(C*CFF,int,idf)
       end if

       if(joper(4).gt.0.and.Msym(is).ne.0) then        ! Z-integrals
        C = C * Msym(is) * (Ssym(is)-1)
        int = Incode_int (7,0,k1,k1,k2,k2)
        Call Iadd_zoef(C,int,idf)
       end if

       End do
       End do
      End do

      End Subroutine ZNO_1ee

!======================================================================
      Subroutine ZNO_2ee (i1,i2,j1,j2)
!======================================================================
!     angular part of matrix elements between two det.w.f.
!     for two-electron operator
!     Calls: Check_boef, Idet_fact, Iadd_zoef.
!----------------------------------------------------------------------
      Use spin_orbitals, only: nsym, kz1,kz2,Lsym,Msym,Ssym, &
                               IPsym1,IPsym2,Isym1,Isym2,nnsym1,nnsym2
      Use BOEF_list,     only: kblk,ncblk, ib_int, boef

      Implicit none
      Integer, intent(in) :: i1,i2,j1,j2
      Integer :: i11,i12,i21,i22,j11,j12,j21,j22
      Integer :: io1,io2,io3,io4, ib1,ib2,ib3,ib4
      Integer :: i,j,k, k1,k2,k3,k4, met,int,idf,kz,ii1,ii2
      Integer :: ibint(4)
      Integer, external :: Idet_fact, Incode_int

!----------------------------------------------------------------------
      if(mod(Lsym(i1)+Lsym(i2)+Lsym(j1)+Lsym(j2),2).ne.0) Return
      if(Msym(i1)+Msym(i2).ne.Msym(j1)+Msym(j2)) Return
      if(Ssym(i1)+Ssym(i2).ne.Ssym(j1)+Ssym(j2)) Return
!----------------------------------------------------------------------

      Call Check_boef(Lsym(i1),Msym(i1),Ssym(i1), &
                      Lsym(i2),Msym(i2),Ssym(i2), &
                      Lsym(j1),Msym(j1),Ssym(j1), &
                      Lsym(j2),Msym(j2),Ssym(j2))

!----------------------------------------------------------------------
!
      i11 = 1;  if(i1.gt.1) i11 = IPsym1(i1-1)+1;  i12 = IPsym1(i1)
      i21 = 1;  if(i2.gt.1) i21 = IPsym1(i2-1)+1;  i22 = IPsym1(i2)
      j11 = 1;  if(j1.gt.1) j11 = IPsym2(j1-1)+1;  j12 = IPsym2(j1)
      j21 = 1;  if(j2.gt.1) j21 = IPsym2(j2-1)+1;  j22 = IPsym2(j2)

      Do k1 = i11,i12;  i=Isym1(k1); ibint(1)=nnsym1(i)
      Do k2 = i21,i22;  j=Isym1(k2); ibint(2)=nnsym1(j)
                        if(k2.le.k1) Cycle

      Do k3 = j11,j12;  i=Isym2(k3);  ibint(3)=nnsym2(i)
      Do k4 = j21,j22;  j=Isym2(k4);  ibint(4)=nnsym2(j)
                        if(k4.le.k3) Cycle

       idf = Idet_fact(k1,k2,k3,k4)

       kz = (-1)**(kz1+kz2+k1+k2+k3+k4)

       ii1 = 1; if(kblk.gt.1) ii1=ncblk(kblk-1)+1; ii2=ncblk(kblk)
       Do i = ii1,ii2
        Call Decode_int (met,k,ib1,ib2,ib3,ib4,IB_int(i))
        io1 = ibint(ib1); io2 = ibint(ib2)
        io3 = ibint(ib3); io4 = ibint(ib4)
        Call Jsym_int(met,io1,io2,io3,io4)
        int = Incode_Int (met,k,io1,io2,io3,io4)
        Call Iadd_zoef(Boef(i)*kz,int,idf)

       End do

      End do;  End do;  End do;  End do

      End Subroutine ZNO_2ee


!======================================================================
      Subroutine Jsym_int(icase,j1,j2,j3,j4)
!======================================================================
!     define "canonical" form of integral
!---------------------------------------------------------------------
      Implicit none
      Integer :: icase,j,j1,j2,j3,j4

      if(icase.eq.3.or.icase.eq.4.or.icase.eq.5) then
       if(j1.gt.j2) then
        j=j1;j1=j2;j2=j; j=j3;j3=j4;j4=j
       elseif(j1.eq.j2.and.j3.gt.j4) then
        j=j3;j3=j4;j4=j
       end if
      end if

      End Subroutine Jsym_int
