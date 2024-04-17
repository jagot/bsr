!=======================================================================
Subroutine DD_OUT
  !=======================================================================
  !     define and output dipole matrix for given initial state
  !-----------------------------------------------------------------------
  Use dbsr_dmat
  Use target_jj, only: etarg
  Use conf_jj

  Implicit none
  Character(64) :: Label1,Label2
  Integer :: i,j,nhm,kch,kpert,ns1,isol,tmp
  Real(8) :: S,SL,SV,E1

  Real(8), allocatable :: Cbra(:,:), Cket(:,:)
  Real(8), allocatable :: BBL(:,:), AAL(:,:), BBV(:,:), AAV(:,:)
  Real(8), allocatable :: CL(:,:),CV(:,:)

  if(ktype.ne.'E') Stop 'dd_out: non-electric-transition case ? '
  if(kpol.ne.1) Stop 'dd_out: kpol <> 1 --> non-dipole case ? '
  !----------------------------------------------------------------------
  !                                      define the initial bound states:
  i=INDEX(BF_b,'.'); AF=BF_b(1:i)//ALS1
  Call Check_file(AF); Open(nub1,file=AF,form='UNFORMATTED')
  rewind(nub1)

  read(nub1) nhm,kch,kpert,ns1,jot1,parity,nstate1

  if(ns1 .ne.ns ) Stop 'dbsr_dmat: ns1 <> ns '
  if(nch1.ne.kch) Stop 'dbsr_dmat: nch1 --> ?'
  if(npert1.ne.kpert) Stop 'dbsr_dmat: npert1 --> ?'
  if(kdm1.ne.nhm) Stop 'dbsr_dmat: nhm1 --> ?'

  ! We make the radial dimension the slow one, since we need to
  ! contract it with the dipole matrix below.
  Allocate(Cket(nstate1,kdm1))
  Do j=1,nstate1
     read(nub1) tmp,Label1
     read(nub1) E1
     read(nub1) Cket(j,:)
  End do

  i=INDEX(BF_b,'.'); AF=BF_b(1:i)//ALS2
  Call Check_file(AF); Open(nub2,file=AF,form='UNFORMATTED')
  rewind(nub2)

  read(nub2) nhm,kch,kpert,ns1,jot2,parity,nstate2

  if(ns1 .ne.ns )       Stop 'dbsr_dmat: ns2 <> ns '
  if(nch2.ne.kch)       Stop 'dbsr_dmat: nch2 --> ?'
  if(npert2.ne.kpert)   Stop 'dbsr_dmat: npert2 --> ?'
  if(kdm2.ne.nhm)       Stop 'dbsr_dmat: nhm2 --> ?'
  if(parity2.ne.parity) Stop 'dbsr_dmat: parity2 --> ?'

  Allocate(Cbra(kdm2,nstate2))
  Do j=1,nstate2
     read(nub2) tmp,Label2
     read(nub2) E1
     read(nub2) Cbra(:,j)
  End do

  ! write(*,*) "Matrix sizes:"
  ! write(*,'(a10,"(",i10,",",i10,")")') "DL", size(DL,1), size(DL,2)
  ! write(*,'(a10,"(",i10,",",i10,")")') "DV", size(DV,1), size(DV,2)
  ! write(*,'(a10,"(",i10,",",i10,")")') "Cbra", size(Cbra,1), size(Cbra,2)
  ! write(*,'(a10,"(",i10,",",i10,")")') "Cket", size(Cket,1), size(Cket,2)

  !----------------------------------------------------------------------
  ! ... calculation and output the dipole matrix:

  ! Cf. bsr/BSR_DMAT3/dd_out.f90
  !-----------------------------------------------------------------------
  ! ... D(kdm1,kdm2) * Cbra(kdm2,nstate2)

  Allocate (BBL(kdm1,nstate2),BBV(kdm1,nstate2))

  BBL = MATMUL (DL,Cbra)
  BBV = MATMUL (DV,Cbra)

  ! ... Cket(nstate1,kdm1) * BB(kdm1,nstate2)

  Allocate (AAL(nstate1,nstate2),AAV(nstate1,nstate2))

  AAL = MATMUL (Cket,BBL)
  AAV = MATMUL (Cket,BBV)

  Deallocate (BBL,BBV)

  Allocate (CL(nstate2,nstate1),CV(nstate2,nstate1))
  CL = 0
  CV = 0

  CL = TRANSPOSE(AAL)
  CV = TRANSPOSE(AAV)

  Deallocate(AAL,AAV)
  !-----------------------------------------------------------------------

  write(AF,'(a,a,a,a)') 'dd.',ALS1,'_',ALS2
  write(*,*) "Dumping dipole matrix between states <", &
       ALS2, "|d|", ALS1, "> to ", AF

  Open(nud,file=AF,form='FORMATTED')
  write(nud, '(a," <",a,"|d|",a,"> ",a)') "! Dipole matrices", ALS2, ALS1, &
       "are dumped linearly, i.e. D(:,1), then D(:,2), and so on."
  write(nud, '(a)') "! First length form, then velocity form."
  write(nud, '(a," (",i10,",",i10,")")') "! Both matrices are of size", nstate2,nstate1
  write(nud, *)
  write(nud, '("ket_sym = ",a)') ALS1
  write(nud, '("bra_sym = ",a)') ALS2
  write(nud, '("nstate1 = ",i10," ! Ket states")') nstate1
  write(nud, '("nstate2 = ",i10," ! Bra states")') nstate2
  write(nud, *)

  write(nud,'("Length form")')
  write(nud,'((5(g24.16)))') CL
  write(nud,'("*")')
  write(nud, *)
  write(nud,'("Velocity form")')
  write(nud,'((5(g24.16)))') CV
  write(nud,'("*")')

  Deallocate(Cbra,Cket,CL,CV)

  close(nud)

End Subroutine dd_out
