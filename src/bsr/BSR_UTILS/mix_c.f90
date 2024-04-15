!======================================================================
!     utility       m i x _ c
!
!                   C O P Y R I G H T -- 2007
!
!     Written by:   Oleg Zatsarinny
!
!======================================================================
!
!     Extract the state expansion from GRASPvu mix-fails
!
!     arguments: 1. name for m-file
!                2. # of seeking solution
!                3. 2J value
!                4. name for result c-file
!                5. eps_c - tolerance for coefficients
!
!     When eps_c > 0, configurations in result c-file are ordered
!     according where weights
!
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)

      Character(6) :: G92MIX
      Character(80) :: AN, AF ,BF
      Character(300) :: AS

      Real(8), allocatable, Dimension(:) :: WT
      Integer, allocatable, Dimension(:) :: IP,JP
      Character(300), allocatable, Dimension(:) :: AC1,AC2,AC3

      Integer(4), Allocatable, Dimension(:) :: IVEC,IATJPO,IASPAR
      Integer(4), Allocatable, Dimension(:) :: ICCMIN

      Integer(4) :: nuc =1       !   name.c
      Integer(4) :: num =2       !   name.m or name.bm
      Integer(4) :: iout=9       !   result.c

!----------------------------------------------------------------------
!                                                           input data:
      iarg = IARGC()
      if(iarg.eq.0) then
       open(iout,file='run_mix_c.bat')
       write(iout,'(a)') 'mix_c 5p5.m 1 0 a.c 0.00000000001'
       write(iout,*)
       write(iout,'(a)') 'rem  You should provide as arguments the following data:'
       write(iout,*)
       write(iout,'(a)') 'rem  name for m-file'
       write(iout,'(a)') 'rem  # of solution'
       write(iout,'(a)') 'rem  [+|-] 2J value'
       write(iout,'(a)') 'rem  name for resulting c-file'
       write(iout,'(a)') 'rem  eps_c - tolerance for coefficients'
       Stop ' '
      elseif(iarg.lt.5) then
       write(*,*) 'You should provide as arguments the following data:'
       write(*,*)
       write(*,*) '1. name for m-file'
       write(*,*) '2. # of solution'
       write(*,*) '3. [+|-] 2J value'
       write(*,*) '4. name for resulting c-file'
       write(*,*) '5. eps_c - tolerance for coefficients'
       Stop ' '
      else
       Call GETARG(1,AN)
       Call GETARG(2,AF); read(AF,*) nn
       Call GETARG(3,AF); read(AF,*) jj
       Call GETARG(4,BF)
       Call GETARG(5,AF); read(AF,*) eps_c
      end if

!----------------------------------------------------------------------
! ... read list of configuration:

      ia=INDEX(AN,'.',BACK=.TRUE.)

      if(AN(ia+1:ia+1).ne.'m') Stop ' mix-file should have extensions m'

      AF=AN(1:ia)//'c'
      Call Check_file(AF)
      Open(nuc,file=AF,STATUS='OLD')

      ncfg=Jdef_ncfg(nuc)
      Allocate(WT(ncfg),AC1(ncfg),AC2(ncfg),AC3(ncfg), &
               JP(ncfg),IP(ncfg))

      rewind(nuc)
      read(nuc,'(a)') AS
      read(nuc,'(a)') AS
      read(nuc,'(a)') AS
      read(nuc,'(a)') AS

      i=0
    1 read(nuc,'(a)',end=2) AS
      if(AS(6:6).ne.'(') go to 1
      i=i+1; AC1(i)=AS
      read(nuc,'(a)') AC2(i)
      read(nuc,'(3x,a)') AC3(i)

      AS=AC3(i); j=LEN_TRIM(AS)
      if(AS(j-2:j-2).eq.'/') then
       read(AS(j-7:j-3),*) JP(i)
      else
       read(AS(j-7:j-1),*) JP(i)
      end if
      JP(i) = JP(i)*2; if(AS(j:j).eq.'-') JP(i)=-JP(i)

      go to 1
    2 Continue

!----------------------------------------------------------------------
! ... read solution:

      Call Check_file(AN)
      Open(num,file=AN,form='UNFORMATTED')
      read(num) G92MIX
      if(G92MIX.ne.'G92MIX') then
       write(*,'(a,a,a)') 'file ',AN,' is not GRASP mix-fail'
       Stop ' '
      end if

!----------------------------------------------------------------------
! ... m-files:

       read(num) NELEC, NCFTOT, NW, NVECTOT, NVECSIZ, NBLOCK

       if(NCFTOT.ne.ncfg) Stop 'ncfg in c-fail <> NCFTOT in mb-fail'
       ic2=0; k = 0
       Do ib = 1,NBLOCK
        read(num) NB, NCFBLK, NEVBLK, IATJP, IASPA
        ic1=ic2+1; ic2=ic2+NCFBLK
        if(jj.ne.(IATJP-1)*IASPA) then
         read(num) i
         read(num) EAV
         read(num) EVEC
        else
         Allocate(ICCMIN(NEVBLK))
         read(num) (ICCMIN(I),I=1,NEVBLK)
         Do i=1,NEVBLK; if(nn.ne.ICCMIN(I)) Cycle; k=i; Exit; End do
         if(k.eq.0) then
          read(num) EAV; read(num) EVEC
         else
          read(num) EAV, (EVAL, I = 1, k)
          EVAL = EVAL + EAV
          read(num) ((WT(i), i=ic1,ic2), J=1,k)
         end if
        end if
        if(k.gt.0) Exit
       End do
       if(k.eq.0) Stop 'Can not find the required solution'

!----------------------------------------------------------------------
! ... order the configurations according their weights:

      Do i=1,ncfg; IP(i)=i; End do

      if(eps_c.gt.0.d0) then
       Do i=ic1,ic2-1
        m=i
        Do j=i+1,ic2
         if(abs(WT(IP(j))).gt.abs(WT(IP(m)))) m=j
        End do
        if(m.ne.i) then
         k=IP(m); IP(m)=IP(i); IP(i)=k
        end if
       End do
      end if

!----------------------------------------------------------------------
! ... output the c-file:

      open(iout,file=BF);  rewind(nuc)

      read(nuc,'(a)') AS
      write(iout,'(a,f16.8)') 'Core subshells:',EVAL

      read(nuc,'(a)') AS;  write(iout,'(a)') TRIM(AS)
      read(nuc,'(a)') AS;  write(iout,'(a)') TRIM(AS)
      read(nuc,'(a)') AS;  write(iout,'(a)') TRIM(AS)
      read(nuc,'(a)') AS;  write(iout,'(a)') TRIM(AS)

      Do ic=ic1,ic2
       i=IP(ic); if(abs(WT(i)).lt.Eps_c) Exit
       AS = AC1(i)
       ii = LEN_TRIM(AS); ii = max(72,ii)
       write(iout,'(a,f12.8)') AS(1:ii),WT(i)
       write(iout,'(a)') AC2(i)
       write(iout,'(3x,a)') AC3(i)
      End do

      write(iout,'(a)') '***'

      End  ! program cfile


!======================================================================
      Integer(4) Function Jdef_ncfg(nu)
!======================================================================
!     defines the number of configuration in c-file (unit nu)
!----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(4), INTENT(in) :: nu
      INTEGER(4) :: ncfg
      CHARACTER(6) :: AS

      rewind(nu)
      ncfg=0
    1 read(nu,'(a)',end=2) AS
      if(AS(1:3).eq.'***') go to 2
      if(AS(6:6).ne.'(') go to 1
      ncfg=ncfg+1
      go to 1
    2 rewind(nu)
      Jdef_ncfg=ncfg

      End Function Jdef_ncfg


!======================================================================
      Subroutine Check_file(AF)
!======================================================================

      Character(*), Intent(in) :: AF

      Logical :: EX

      Inquire(file=AF,exist=EX)

      if(.not.EX) then

       write(*,*) ' can not find file  ',AF;  Stop

      end if

      End Subroutine Check_file


