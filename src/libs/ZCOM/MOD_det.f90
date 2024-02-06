!----------------------------------------------------------------------
!     Subroutines for processing overlap determinants:
!
!     DET_list,  alloc_det, Iadd_det
!     DEF_list,  alloc_def, Iadd_def
!
!     NDET_list,  alloc_ndet, Iadd_ndet
!     NDEF_list,  alloc_ndef, Iadd_ndef
!
!     NDET_IDET, IDET_SIMP
!
!     DET_list and DEF_list contain the overlap determinants
!     and overlap factors for all configuration symmetries,
!     whereas NDET_list and NDEF_list contain the information
!     only for two configuration under consideration.
!     These additional NDET, NDEF lists were introduced in hope
!     to reduce the seeking time in large common lists DET and DEF.
!
!     All four modules DET_list, DEF_list, NDET_list, NDEF_list
!     have the identical structure and  differ only by names
!     for variables.
!
!     The connection between DET,DEF and NDET,NDEF lists is given
!     by subroutine NDET_IDET.
!
!     Each overlap determinant is codded as list of involved orbitals:
!     i1*base+j1, i2*base+j2, ... , where {i},{j} - pointers to orbitals,
!     in left and right bra-vectors, respectively.
!
!----------------------------------------------------------------------


!======================================================================
      Module DET_list
!======================================================================
!     Containes the overlap determinants for all symmetries
!----------------------------------------------------------------------

      Implicit none

      Integer :: ndet = 0       ! number of determinants
      Integer :: mdet = 0       ! current dimension for KPD,IPD
      Integer :: idet = 2**15   ! initial dimension
      Integer :: jdet = 2**2    ! average size of one determinant
      Integer :: kdet = 0       ! current dimension for NPD
      Integer :: ldet = 0       ! last filled element in NPD
      Integer :: jmdet = 2**2   ! maximim size of one determinant


      Integer :: ibd = 2**15    ! overlap determinant packing basis

      Integer, allocatable :: KPD(:),IPD(:),NPD(:),JPD(:)

      ! KPD(i) - dimension of i-th determinant
      ! IPD(i) - its pointer in the list NPD
      ! JPD(i) - ordering pointer

      End Module DET_list


!======================================================================
      Subroutine alloc_det(m)
!======================================================================
      Use det_list

      Implicit none
      Integer, Intent(in) :: m
      Integer, Allocatable :: ARR(:)

      if(m.le.0) then
       if(allocated(KPD)) Deallocate (KPD,IPD,NPD,JPD)
       mdet=0; kdet=0; ldet=0; ndet=0
       if(m.lt.0) then
        mdet=idet; kdet = mdet*jdet
        Allocate(KPD(mdet),IPD(mdet),NPD(kdet),JPD(mdet))
       end if
      elseif(.not.allocated(KPD)) then
       mdet=m; kdet=mdet*jdet
       Allocate(KPD(mdet),IPD(mdet),NPD(kdet),JPD(mdet))
      elseif(m.le.mdet) then
       Return
      elseif(ndet.eq.0) then
       Deallocate (KPD,IPD,NPD,JPD)
       mdet=m; kdet=mdet*jdet; ldet=0
       Allocate(KPD(mdet),IPD(mdet),NPD(kdet),JPD(mdet))
      else
       jdet=ldet/ndet+1; mdet=m; kdet=mdet*jdet
       if(kdet.lt.ldet+10*jdet) kdet=ldet+10*jdet
       Allocate(ARR(ldet))
       ARR(1:ndet)=KPD(1:ndet); Deallocate(KPD)
       Allocate(KPD(mdet)); KPD(1:ndet)=ARR(1:ndet)
       ARR(1:ndet)=IPD(1:ndet); Deallocate(IPD)
       Allocate(IPD(mdet)); IPD(1:ndet)=ARR(1:ndet)
       ARR(1:ldet)=NPD(1:ldet); Deallocate(NPD)
       Allocate(NPD(kdet)); NPD(1:ldet)=ARR(1:ldet)
       ARR(1:ndet)=JPD(1:ndet); Deallocate(JPD)
       Allocate(JPD(mdet)); JPD(1:ndet)=ARR(1:ndet)
       Deallocate(ARR)
       write(*,*) 'realloc_det : m = ',m,idet
      end if

      End Subroutine alloc_DET


!======================================================================
      Integer Function Iadd_det  (kd,NP)
!======================================================================
!     add new overlap determinant to DET_list
!----------------------------------------------------------------------
      Use det_list

      Implicit none
      Integer , Intent(in) :: kd, NP(kd)
      Integer :: i,j,k,m,ip,i1,i2

      Iadd_det  = 0
      if(kd.le.0) Return
      if(mdet.eq.0) Call Alloc_DET (idet)
      if(kd.gt.jmdet) jmdet = kd

! ... check if the same det. is already in the list:

      i1=1; i2=ndet
    1 if(i1.gt.i2) go to 2
      i=(i1+i2)/2; j=jpd(i)
      if    (kd.lt.kpd(j)) then;  i2 = i - 1
      elseif(kd.gt.kpd(j)) then;  i1 = i + 1
      else
       ip = IPD(j); m = 0
       Do k = 1,kd; ip=ip+1
        if(NP(k).eq.NPD(ip)) Cycle
        if(NP(k).lt.NPD(ip)) then; m = -1; Exit; end if
        if(NP(k).gt.NPD(ip)) then; m = +1; Exit; end if
       End do
       if(m.eq.0) then
        Iadd_det  = j;  Return      ! or i ???
       elseif(m.eq.-1) then
        i2 = i - 1
       else
        i1 = i + 1
       end if
      end if
      go to 1
    2 Continue

! ... Add new determinant:

      if(ndet.ge.mdet.or.ldet+kd.gt.kdet) Call Alloc_det (mdet+idet)
      ndet=ndet+1
      KPD(ndet)=kd
      IPD(ndet)=ldet
      Do i=1,kd; ldet=ldet+1; NPD(ldet)=NP(i); End do

      Do i=ndet,i1+1,-1; jpd(i)=jpd(i-1); End do
      jpd(i1)=ndet
      Iadd_det =ndet     ! or ndet ???

      End Function Iadd_det


!======================================================================
      Subroutine Record_det (nu)
!======================================================================
!     Containes the overlap determinants for all symmetries
!----------------------------------------------------------------------
      Use det_list

      Implicit none

      Integer :: nu, i

      write(nu) ndet,ldet,jmdet
      if(ndet.eq.0) Return

      Do i = 1,ndet
       write(nu) kpd(i),ipd(i),jpd(i),NPD(ipd(i)+1:ipd(i)+kpd(i))
      End do

      End Subroutine Record_det


!======================================================================
      Subroutine Read_det (nu)
!======================================================================
!     Containes the overlap determinants for all symmetries
!----------------------------------------------------------------------
      Use det_list

      Implicit none
      Integer :: nu, i

      read(nu) ndet,ldet,jmdet
      if(ndet.eq.0) Return

      if(allocated(KPD)) Deallocate (KPD,IPD,NPD,JPD)
      mdet = ndet; kdet = ldet
      Allocate(KPD(mdet),IPD(mdet),NPD(kdet),JPD(mdet))

      Do i = 1,ndet
       read(nu) kpd(i),ipd(i),jpd(i),NPD(ipd(i)+1:ipd(i)+kpd(i))
      End do

      End Subroutine Read_det

!======================================================================
      Subroutine Load_det (nu)
!======================================================================
!     Containes the overlap determinants for all symmetries
!----------------------------------------------------------------------
      Use det_list

      Implicit none
      Integer :: nu

      read(nu) ndet,ldet,jmdet
      if(ndet.eq.0) Return

      if(allocated(kpd)) Deallocate (KPD,IPD,NPD,JPD)
      mdet = ndet; kdet = ldet
      Allocate(KPD(mdet),IPD(mdet),NPD(kdet),JPD(mdet))

      read(nu) kpd
      read(nu) ipd
      read(nu) npd
      read(nu) jpd

      End Subroutine Load_det



!======================================================================
      Subroutine Write_det (nu)
!======================================================================
!     Containes the overlap determinants for all symmetries
!----------------------------------------------------------------------
      Use det_list

      Implicit none
      Integer :: nu

      write(nu) ndet,ldet,jmdet
      if(ndet.eq.0) Return

      write(nu) kpd(1:ndet)
      write(nu) ipd(1:ndet)
      write(nu) npd(1:ldet)
      write(nu) jpd(1:ndet)

      End Subroutine Write_det



!======================================================================
      Module DEF_list
!======================================================================
!     Containes the overlap factors, as the list of the number of
!     involved overlap determinants and their positions in the
!     common det_list.
!
!     KPF(i) - number of det.s in i-th overlap factor
!     IPF(i) - pointer on the list of corr. det.s in the NPF
!     NPF(ip+1:ip+kd) - list of pointers on det.s in the DET_list
!                       and their powers
!     JPD(i) - ordering pointer
!----------------------------------------------------------------------
      Implicit none

      Integer :: ndef = 0       ! number of determinants
      Integer :: mdef = 0       ! current dimentsion of list
      Integer :: idef = 2**15   ! supposed max. dimentsion
      Integer :: jdef = 2**3    ! average number of det.s
      Integer :: kdef = 0       ! dimension of all def.s
      Integer :: ldef = 0       ! dimension of all def.s
      Integer :: jmdef = 2*3    ! maximim size of one determinant

      Integer :: ibf = 2**4     ! overlap factors packing basis

      Integer, Allocatable :: KPF(:),IPF(:),NPF(:),JPF(:)

      End Module DEF_list


!======================================================================
      Subroutine alloc_def (m)
!======================================================================
      Use def_list

      Implicit none
      Integer, Intent(in) :: m
      Integer, Allocatable :: ARR(:)

      if(m.le.0) then
       if(allocated(KPF)) Deallocate(KPF,IPF,NPF,JPF)
       mdef = 0; kdef = 0; ndef =0; ldef = 0
       if(m.lt.0) then
        mdef=idef; kdef = mdef*jdef
        Allocate(KPF(mdef),IPF(mdef),NPF(kdef),JPF(mdef))
       end if
      elseif(.not.allocated(KPF)) then
       mdef = m; kdef = mdef*jdef
       Allocate(KPF(mdef),IPF(mdef),NPF(kdef),JPF(mdef))
      elseif(m.le.mdef) then
       Return
      elseif(ndef.eq.0) then
       Deallocate (KPF,IPF,NPF,JPF)
       mdef = m; kdef = mdef*jdef
       Allocate(KPF(mdef),IPF(mdef),NPF(kdef),JPF(mdef))
      else
       jdef=ldef/ndef+1; mdef=m; kdef=mdef*jdef
       if(kdef.lt.ldef+10*jdef) kdef=ldef+10*jdef
       Allocate(ARR(ldef))
       ARR(1:ndef)=KPF(1:ndef); Deallocate(KPF)
       Allocate(KPF(mdef)); KPF(1:ndef)=ARR(1:ndef)
       ARR(1:ndef)=IPF(1:ndef); Deallocate(IPF)
       Allocate(IPF(mdef)); IPF(1:ndef)=ARR(1:ndef)
       ARR(1:ndef)=JPF(1:ndef); Deallocate(JPF)
       Allocate(JPF(mdef)); JPF(1:ndef)=ARR(1:ndef)
       ARR(1:ldef)=NPF(1:ldef); Deallocate(NPF)
       Allocate(NPF(kdef)); NPF(1:ldef)=ARR(1:ldef)
       Deallocate(ARR)
       write(*,*) 'realloc_def : m = ',m,idef
      end if

      End Subroutine alloc_def


!======================================================================
      Integer Function Iadd_def  (kd,NP)
!======================================================================
!     add new overlap determinant to DET_list
!----------------------------------------------------------------------
      Use def_list

      Implicit none
      Integer , Intent(in) :: kd, NP(kd)
      Integer :: i,j,k,m,ip,i1,i2

      Iadd_def  = 0
      if(kd.le.0) Return
      if(mdef.eq.0) Call Alloc_DEF (idef)
      if(kd.gt.jmdef) jmdef=kd

! ... check if the same det. is already in the list:

      i1=1; i2=ndef
    1 if(i1.gt.i2) go to 2
      i=(i1+i2)/2; j=jpf(i)
      if    (kd.lt.kpf(j)) then;  i2 = i - 1
      elseif(kd.gt.kpf(j)) then;  i1 = i + 1
      else
       ip = IPF(j); m = 0
       Do k = 1,kd; ip=ip+1
        if(NP(k).eq.NPF(ip)) Cycle
        if(NP(k).lt.NPF(ip)) then; m = -1; Exit; end if
        if(NP(k).gt.NPF(ip)) then; m = +1; Exit; end if
       End do
       if(m.eq.0) then
        Iadd_def  = j; Return      ! or i ???
       elseif(m.eq.-1) then
        i2 = i - 1
       else
        i1 = i + 1
       end if
      end if
      go to 1
    2 Continue

! ... Add new determinant:

      if(ndef.ge.mdef.or.ldef+kd.gt.kdef) Call Alloc_DEF (mdef+idef)
      ndef=ndef+1
      KPF(ndef)=kd
      IPF(ndef)=ldef
      Do i=1,kd; ldef=ldef+1; NPF(ldef)=NP(i); End do

      Do i=ndef,i1+1,-1; jpf(i)=jpf(i-1); End do
      jpf(i1)=ndef
      Iadd_def =ndef     ! or i1 ???

      End Function Iadd_def


!======================================================================
      Subroutine Record_def (nu)
!======================================================================
!     Containes the overlap determinants for all symmetries
!----------------------------------------------------------------------
      Use def_list

      Implicit none
      Integer :: nu, i

      write(nu) ndef,ldef,jmdef
      if(ndef.eq.0) Return

      Do i = 1,ndef
       write(nu) kpf(i),ipf(i),jpf(i),npf(ipf(i)+1:ipf(i)+kpf(i))
      End do

      End Subroutine Record_def


!======================================================================
      Subroutine Write_def (nu)
!======================================================================
!     Containes the overlap determinants for all symmetries
!----------------------------------------------------------------------
      Use def_list

      Implicit none
      Integer :: nu

      write(nu) ndef,ldef,jmdef
      if(ndef.eq.0) Return

      write(nu) kpf(1:ndef)
      write(nu) ipf(1:ndef)
      write(nu) npf(1:ldef)
      write(nu) jpf(1:ndef)

      End Subroutine Write_def


!======================================================================
      Subroutine Load_def (nu)
!======================================================================
!     Containes the overlap determinants for all symmetries
!----------------------------------------------------------------------
      Use def_list

      Implicit none
      Integer :: nu

      read(nu) ndef,ldef,jmdef
      if(ndef.eq.0) Return

      if(allocated(kpf)) Deallocate (KPF,IPF,NPF,JPF)
      mdef = ndef; kdef = ldef
      Allocate(KPF(mdef),IPF(mdef),NPF(kdef),JPF(mdef))

      read(nu) kpf
      read(nu) ipf
      read(nu) npf
      read(nu) jpf

      End Subroutine Load_def


!======================================================================
      Subroutine Read_def (nu)
!======================================================================
!     Containes the overlap determinants for all symmetries
!----------------------------------------------------------------------
      Use def_list

      Implicit none

      Integer :: nu, i

      read(nu) ndef,ldef,jmdef
      if(ndef.eq.0) Return

      if(allocated(kpf)) Deallocate (KPF,IPF,NPF,JPF)
      mdef = ndef; kdef = ldef
      Allocate(KPF(mdef),IPF(mdef),NPF(kdef),JPF(mdef))

      Do i = 1,ndef
       read(nu) kpf(i),ipf(i),jpf(i),npf(ipf(i)+1:ipf(i)+kpf(i))
      End do

      End Subroutine Read_def


!======================================================================
      Module NDET_list
!======================================================================
!     Containes the overlap determinants for all symmetries
!----------------------------------------------------------------------
      Implicit none

      Integer :: ndet = 0       ! number of determinants
      Integer :: mdet = 0       ! current dimension for KPD,IPD
      Integer :: idet = 2**15   ! initial dimension
      Integer :: jdet = 2**2    ! average size of one determinant
      Integer :: kdet = 0       ! current dimension for NPD
      Integer :: ldet = 0       ! last filled element in NPD

      Integer :: ibd = 2**15    ! overlap determinant packing basis

      Integer, allocatable :: KPD(:),IPD(:),NPD(:),JPD(:)

      ! KPD(i) - dimension of i-th determinant
      ! IPD(i) - its pointer in the list NPD
      ! JPD(i) - ordering pointer

      End Module NDET_list


!======================================================================
      Subroutine alloc_ndet(m)
!======================================================================
      Use ndet_list

      Implicit none
      Integer, Intent(in) :: m
      Integer, Allocatable :: ARR(:)

      if(m.le.0) then
       if(allocated(KPD)) Deallocate (KPD,IPD,NPD,JPD)
       mdet=0; kdet=0; ldet=0; ndet=0
       if(m.lt.0) then
        mdet=idet; kdet = mdet*jdet
        Allocate(KPD(mdet),IPD(mdet),NPD(kdet),JPD(mdet))
       end if
      elseif(.not.allocated(KPD)) then
       mdet=m; kdet=mdet*jdet
       Allocate(KPD(mdet),IPD(mdet),NPD(kdet),JPD(mdet))
      elseif(m.le.mdet) then
       Return
      elseif(ndet.eq.0) then
       Deallocate (KPD,IPD,NPD,JPD)
       mdet=m; kdet=mdet*jdet; ldet=0
       Allocate(KPD(mdet),IPD(mdet),NPD(kdet),JPD(mdet))
      else
       jdet=ldet/ndet+1; mdet=m; kdet=mdet*jdet
       if(kdet.lt.ldet+10*jdet) kdet=ldet+10*jdet
       Allocate(ARR(ldet))
       ARR(1:ndet)=KPD(1:ndet); Deallocate(KPD)
       Allocate(KPD(mdet)); KPD(1:ndet)=ARR(1:ndet)
       ARR(1:ndet)=IPD(1:ndet); Deallocate(IPD)
       Allocate(IPD(mdet)); IPD(1:ndet)=ARR(1:ndet)
       ARR(1:ldet)=NPD(1:ldet); Deallocate(NPD)
       Allocate(NPD(kdet)); NPD(1:ldet)=ARR(1:ldet)
       ARR(1:ndet)=JPD(1:ndet); Deallocate(JPD)
       Allocate(JPD(mdet)); JPD(1:ndet)=ARR(1:ndet)
       Deallocate(ARR)
       write(*,*) 'realloc_ndet : m = ',m,idet
      end if

      End Subroutine alloc_NDET


!======================================================================
      Integer Function Iadd_ndet  (kd,NP)
!======================================================================
!     add new overlap determinant to NDET_list
!----------------------------------------------------------------------
      Use ndet_list

      Implicit none
      Integer , Intent(in) :: kd, NP(kd)
      Integer :: i,j,k,m,ip,i1,i2

      Iadd_ndet  = 0
      if(kd.le.0) Return
      if(mdet.eq.0) Call Alloc_NDET (idet)

! ... check if the same det. is already in the list:

      i1=1; i2=ndet
    1 if(i1.gt.i2) go to 2
      i=(i1+i2)/2; j=jpd(i)
      if    (kd.lt.kpd(j)) then;  i2 = i - 1
      elseif(kd.gt.kpd(j)) then;  i1 = i + 1
      else
       ip = IPD(j); m = 0
       Do k = 1,kd; ip=ip+1
        if(NP(k).eq.NPD(ip)) Cycle
        if(NP(k).lt.NPD(ip)) then; m = -1; Exit; end if
        if(NP(k).gt.NPD(ip)) then; m = +1; Exit; end if
       End do
       if(m.eq.0) then
        Iadd_ndet  = j;  Return      ! or i ???
       elseif(m.eq.-1) then
        i2 = i - 1
       else
        i1 = i + 1
       end if
      end if
      go to 1
    2 Continue

! ... Add new determinant:

      if(ndet.ge.mdet.or.ldet+kd.gt.kdet) Call Alloc_ndet (mdet+idet)
      ndet=ndet+1
      KPD(ndet)=kd
      IPD(ndet)=ldet
      Do i=1,kd; ldet=ldet+1; NPD(ldet)=NP(i); End do

      Do i=ndet,i1+1,-1; jpd(i)=jpd(i-1); End do
      jpd(i1)=ndet
      Iadd_ndet =ndet     ! or ndet ???

      End Function Iadd_ndet


!======================================================================
      Module NDEF_list
!======================================================================
!     Containes the overlap factors, as the list of the number of
!     involved overlap determinants and their positions in the
!     common det_list.
!
!     KPF(i) - number of det.s in i-th overlap factor
!     IPF(i) - pointer on the list of corr. det.s in the NPF
!     NPF(ip+1:ip+kd) - list of pointers on det.s in the DET_list
!                       and their powers
!     JPD(i) - ordering pointer
!----------------------------------------------------------------------
      Implicit none

      Integer :: ndef = 0       ! number of determinants
      Integer :: mdef = 0       ! current dimentsion of list
      Integer :: idef = 2**15   ! supposed max. dimentsion
      Integer :: jdef = 2**3    ! average number of det.s
      Integer :: kdef = 0       ! dimension of all def.s
      Integer :: ldef = 0       ! dimension of all def.s

      Integer :: ibf = 2**4     ! overlap factors packing basis

      Integer, Allocatable :: KPF(:),IPF(:),NPF(:),JPF(:)

      End MODULE NDEF_list


!======================================================================
      Subroutine alloc_ndef (m)
!======================================================================
      Use ndef_list

      Implicit none
      Integer, Intent(in) :: m
      Integer, Allocatable :: ARR(:)

      if(m.le.0) then
       if(allocated(KPF)) Deallocate(KPF,IPF,NPF,JPF)
       mdef = 0; kdef = 0; ndef =0; ldef = 0
       if(m.lt.0) then
        mdef=idef; kdef = mdef*jdef
        Allocate(KPF(mdef),IPF(mdef),NPF(kdef),JPF(mdef))
       end if
      elseif(.not.allocated(KPF)) then
       mdef = m; kdef = mdef*jdef
       Allocate(KPF(mdef),IPF(mdef),NPF(kdef),JPF(mdef))
      elseif(m.le.mdef) then
       Return
      elseif(ndef.eq.0) then
       Deallocate (KPF,IPF,NPF,JPF)
       mdef = m; kdef = mdef*jdef
       Allocate(KPF(mdef),IPF(mdef),NPF(kdef),JPF(mdef))
      else
       jdef=ldef/ndef+1; mdef=m; kdef=mdef*jdef
       if(kdef.lt.ldef+10*jdef) kdef=ldef+10*jdef
       Allocate(ARR(ldef))
       ARR(1:ndef)=KPF(1:ndef); Deallocate(KPF)
       Allocate(KPF(mdef)); KPF(1:ndef)=ARR(1:ndef)
       ARR(1:ndef)=IPF(1:ndef); Deallocate(IPF)
       Allocate(IPF(mdef)); IPF(1:ndef)=ARR(1:ndef)
       ARR(1:ndef)=JPF(1:ndef); Deallocate(JPF)
       Allocate(JPF(mdef)); JPF(1:ndef)=ARR(1:ndef)
       ARR(1:ldef)=NPF(1:ldef); Deallocate(NPF)
       Allocate(NPF(kdef)); NPF(1:ldef)=ARR(1:ldef)
       Deallocate(ARR)
       write(*,*) 'realloc_ndef : m = ',m,idef
      end if

      End Subroutine alloc_ndef


!======================================================================
      Integer Function Iadd_ndef  (kd,NP)
!======================================================================
!     add new overlap determinant to DET_list
!----------------------------------------------------------------------
      Use ndef_list

      Implicit none
      Integer , Intent(in) :: kd, NP(kd)
      Integer :: i,j,k,m,ip,i1,i2

      Iadd_ndef  = 0
      if(kd.le.0) Return
      if(mdef.eq.0) Call Alloc_NDEF (idef)

! ... check if the same det. is already in the list:

      i1=1; i2=ndef
    1 if(i1.gt.i2) go to 2
      i=(i1+i2)/2; j=jpf(i)
      if    (kd.lt.kpf(j)) then;  i2 = i - 1
      elseif(kd.gt.kpf(j)) then;  i1 = i + 1
      else
       ip = IPF(j); m = 0
       Do k = 1,kd; ip=ip+1
        if(NP(k).eq.NPF(ip)) Cycle
        if(NP(k).lt.NPF(ip)) then; m = -1; Exit; end if
        if(NP(k).gt.NPF(ip)) then; m = +1; Exit; end if
       End do
       if(m.eq.0) then
        Iadd_ndef  = j;  Return      ! or i ???
       elseif(m.eq.-1) then
        i2 = i - 1
       else
        i1 = i + 1
       end if
      end if
      go to 1
    2 Continue

! ... Add new determinant:

      if(ndef.ge.mdef.or.ldef+kd.gt.kdef) Call Alloc_NDEF (mdef+idef)
      ndef=ndef+1
      KPF(ndef)=kd
      IPF(ndef)=ldef
      Do i=1,kd; ldef=ldef+1; NPF(ldef)=NP(i); End do

      Do i=ndef,i1+1,-1; jpf(i)=jpf(i-1); End do
      jpf(i1)=ndef
      Iadd_ndef =ndef

      End Function Iadd_ndef


!======================================================================
      Subroutine Ndet_Idet
!======================================================================
!     Add NDET and NDEF lists to the DET and DEF lists.
!     Connection is given in IPF array.
!     The NDET and NDEF lists are then nulified.
!----------------------------------------------------------------------
      Use ndet_list
      Use ndef_list

      Implicit none
      Integer :: i,j,ip,jp,id,kd,ns
      Integer, External :: Iadd_det, Iadd_def, ISORT

      if(ndet.le.0) Return
      Do id=1,ndet
       kd=KPD(id); ip=IPD(id); IPD(id)=Iadd_det(kd,NPD(ip+1))
      End do
      ndet = 0; ldet = 0

      if(ndef.le.0) Return
      Do id=1,ndef
       kd=KPF(id); ip=IPF(id)
       Do i=ip+1,ip+kd
        j=NPF(i)/ibf; ns=mod(NPF(i),ibf); jp=IPD(j)
        NPF(i) = jp*ibf + ns
       End do
       i = ISORT (kd,NPF(ip+1))
       IPF(id) = Iadd_def (kd,NPF(ip+1))
      End do
      ndef = 0; ldef = 0

      End Subroutine Ndet_Idet


!======================================================================
      Integer Function IDET_SIMP(kz,kn,N1,N2)
!======================================================================
!     simplify the determinant kn,N1,N2
!     accoding the orthogonality conditions for radial w.f.'s
!     kz - number of needed permutations
!     IDET_SIMP = 0,1,2 with overlap determinant = 0,1 or some value
!----------------------------------------------------------------------
      Implicit none
      Integer, Intent(inout) :: kz,kn
      Integer, Intent(inout) :: N1(*),N2(*)
      Integer :: i,ii,i1,i2, k,kk,k1,k2, m1,m2
      Integer, external :: IORT

      if(kn.le.0) Stop ' IDET_SIMP: kn <= 0'

      IDET_SIMP=0
!----------------------------------------------------------------------
!                       Check for a row with only one non-zero element:
    1  Do i1=1,kn
       k=0
       Do i2=1,kn
        m1=max(N1(i1),N2(i2)); m2=min(N1(i1),N2(i2)); ii=IORT(m1,m2)
        if(ii.ne.0) then
         k=k+1; kk=ii;  if(k.gt.1.or.kk.ne.1) Exit
         k1=i1; k2=i2
        end if
       End do
       if(k.eq.0) Return; if(k.eq.1.and.kk.eq.1) go to 2
      End do

!----------------------------------------------------------------------
!                   Check for a colum with only one non-zero element:

      Do i2=1,kn
       k=0
       Do i1=1,kn
         m1=max(N1(i1),N2(i2)); m2=min(N1(i1),N2(i2)); ii=IORT(m1,m2)
        if(ii.ne.0) then
         k=k+1; kk=ii; if(k.gt.1.or.kk.ne.1) Exit
         k1=i1; k2=i2
        end if
       End do
       if(k.eq.0) Return; if(k.eq.1.and.kk.eq.1) go to 2
      End do

      go to 3
!-----------------------------------------------------------------------
!                                                 the case of <k1|k2>=1:
    2 kn=kn-1
      if(kn.eq.0) then
       IDET_SIMP=1; Return
      end if
      kz=kz+k1+k2
      Do i=k1,kn; N1(i)=N1(i+1); End do
      Do i=k2,kn; N2(i)=N2(i+1); End do
      go to 1
!-----------------------------------------------------------------------
!                                                  ordering of elements:
    3 Continue

      Do i1=1,kn-1
       Do i2=i1+1,kn
        if(N1(i1).gt.N1(i2)) then
         kk=N1(i1);  N1(i1)=N1(i2); N1(i2)=kk; kz=kz+1
        end if
        if(N2(i1).gt.N2(i2)) then
         kk=N2(i1);  N2(i1)=N2(i2); N2(i2)=kk; kz=kz+1
        end if
       End do
      End do

      IDET_SIMP=2

      End Function IDET_SIMP
