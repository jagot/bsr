!======================================================================
!     utility       C F I L E
!
!                   C O P Y R I G H T -- 2004
!
!     Written by:   Oleg Zatsarinny
!
!======================================================================
!
!     Extract the state from  l(j)-file in separate c-files
!
!     arguments: 1. name for l- or j-file
!                2. index of solution as in l(j)-file
!                3. 2J value (for j-files)
!                4. name for result c-file
!                5. eps_c - tolerance for weights
!
!     When eps_c > 0, configurations in result c-file are ordered
!     according to their weights
!
!----------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)

      Logical :: EX
      Character(80) :: AS, clousd
      Character(40) :: AN, AF ,BF

      Real(8), Allocatable :: WT(:)
      Integer, Allocatable :: IP(:)
      Character(64), Allocatable :: AC(:)
      Character(60), Allocatable :: AU(:)

      Integer :: nuc =1       !   name.c
      Integer :: nuj =2       !   name.l or name.j
      Integer :: iout=3       !   result.c

!----------------------------------------------------------------------
!                                                           input data:
      iarg = IARGC()

      if(iarg.lt.5) then
        write(*,'(/a)') 'cfile extracts the solution from l(j)-file in separate c-file'
        write(*,'(/a)') 'Call as:    cfile name.l(j) index 2J name.c eps_c'
        write(*,'(/a)') 'name.l(j) - input file with solutions'
        write(*,'( a)') 'index     - index of the solution in l(j)-file '
        write(*,'( a)') '2J        - 2J-value in j-files, 0 otherwise'
        write(*,'( a)') 'name      - name for output c-file and w-files'
        write(*,'( a)') 'eps_c     - tolerance for weights'
        write(*,'(/a)') 'When eps_c > 0, configurations in result c-file are ordered'
        write(*,'( a)') 'according to their weights'
        Stop ' '
      else
        Call GETARG(1,AN)
        Call GETARG(2,AF); read(AF,*) nn
        Call GETARG(3,AF); read(AF,*) jj
        Call GETARG(4,BF)
        Call GETARG(5,AF); read(AF,*) eps_c
      end if

      i=LEN_TRIM(BF); if(i.gt.2.and.BF(i-1:i).ne.'.c')  BF=BF(1:I)//'.c'

!----------------------------------------------------------------------

      i=LEN_TRIM(AN); if(AN(i:i).eq.'l') jj=-1

! ... read list of configuration:

      AF=AN;  AF(i:i)='c'
      Open(nuc,file=AF,STATUS='OLD')

      ncfg=Idef_ncfg(nuc)
      if(allocated(WT)) Deallocate(WT,AC,AU,IP)
      Allocate(WT(ncfg), AC(ncfg), AU(ncfg), IP(ncfg))

      rewind(nuc)
      read(nuc,'(a)') AS
      read(nuc,'(a)') Clousd

      i=0
    2 read(nuc,'(a)',end=3) AS
      if(AS(1:1).eq.'*') go to 3
      if(AS(5:5).ne.'(') go to 2
      i=i+1; AC(i)=AS(1:64)
      read(nuc,'(a)') AU(i)
      go to 2
    3 Close(nuc)

! ... read solution:

      Open(nuj,file=AN,STATUS='OLD')

      rewind(nuj)
    5 read(nuj,'(a)',end=9) AS
      if(AS(3:5).ne.'2*J') go to 5
      read(AS,'(7x,i5,10x,i4)') j,nj

      if(jj.ge.0.and.j.ne.jj) go to 5
      Do is=1,nj
       read(nuj,'(/i6, f16.8,i6)') n,E
       read(nuj,'(7F11.8)') (WT(i),i=1,ncfg)
       if(n.eq.nn) Exit
      End do

! ... order the configurations according their weights:

      Do i=1,ncfg; IP(i)=i; End do

      if(eps_c.gt.0.d0) then
      Do i=1,ncfg-1
       m=i
       Do j=i+1,ncfg
        if(abs(WT(IP(j))).gt.abs(WT(IP(m)))) m=j
       End do
       if(m.ne.i) then
        k=IP(m); IP(m)=IP(i); IP(i)=k
       end if
      End do
      end if

! ... output the c-file:

      open(iout,file=BF)
      if(jj.ge.0) then
       write(iout,'(15x,f16.8,a,i3)') E,'   2*J = ',jj
      else
       write(iout,'(15x,f16.8)') E
      end if
      write(iout,'(a)') Clousd

      Do m=1,ncfg
       i=IP(m); if(abs(WT(i)).lt.Eps_c) Exit
       write(iout,'(a64,f11.8)') AC(i),WT(i)
       AS=AU(i); ii=len_trim(AS)
       write(iout,'(a)') AS(1:ii)
      End do
      write(iout,'(a)') '*'
      Close(iout)

       i = LEN_TRIM(AN); AN(i:i)='w'
       if(Icheck_file(AN).eq.1) then
        i = LEN_TRIM(BF); BF(i:i)='w'
        write(AS,'(a,a,a,a)') 'cp ',trim(AN),' ',trim(BF)
        i = SYSTEM(AS)
       end if

      go to 10

    9 write(*,*) ' No solution for jj,n =', jj,nn
   10 Continue

      End  ! program cfile

!======================================================================
      Integer Function Idef_ncfg(nu)
!======================================================================
!
!     gives the number of configuration in c-file (unit nu)
!
!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(in) :: nu

      INTEGER :: ncfg
      Character(5) :: AS

      rewind(nu)
      ncfg=0
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      read(nu,'(a)') AS
      ncfg=ncfg+1
      go to 1
    2 rewind(nu)
      Idef_ncfg=ncfg

      End Function Idef_ncfg



!======================================================================
      Integer Function Icheck_file(AF)
!======================================================================
!     check if the file AF exists
!----------------------------------------------------------------------
      Character(*), Intent(in) :: AF
      Logical :: EX
      Inquire(file=AF,exist=EX)
      Icheck_file = 1
      if(.not.EX) Icheck_file = 0
      End Function Icheck_file

