!======================================================================
!
!     utility     Z G E N C O N F
!
!                 C O P Y R I G H T -- 2004
!
!     Written by:   Oleg Zatsarinny
!
!======================================================================
!
!     preperation of list of possible configurations from list of
!     electron ocupations:
!
!     electron.inp  -->  conf.inp
!
!     see electron.inp for more details
!
!----------------------------------------------------------------------

      Use conf_LS; Use orb_LS

      Character(80) :: AS=' '
      
      Integer, Allocatable :: iq_min(:),iq_max(:)

      Integer :: nu1=1; Character(40) :: AF_inp = 'electron.inp'
      Integer :: nu2=2; Character(40) :: AF_out = 'conf.inp'

!---------------------------------------------------------------------

      mo=0; Call read_iarg('mo',mo); if(mo.eq.0) mo=msh
      Call read_aarg('inp',AF_inp)
      Call read_aarg('out',AF_out)

      iarg = command_argument_count() 
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AS)
      if(AS.eq.'?'.or.Icheck_file(AF_inp).eq.0) then
        write(*,'(/a)') 'zgenconf generates list of configurations (conf.inp)'
        write(*,'(/a)') 'based on the allowed occupation numbers (electron.inp)'
        write(*,'(/a)') 'Call as:  zgenconf [mo=..]'
        write(*,'(/a)') 'Input:    electron.inp (will be created if absent)'
        write(*,'(/a)') 'Results:  conf.inp'
        write(*,'(/a)') 'mo - optional restriction on the number of shells'
        if(Icheck_file(AF_inp).eq.0) Call Create_electron_inp 
        write(*,'(/a)') 
        Stop
      end if

!---------------------------------------------------------------------

      open(nu1,file=AF_inp,status='OLD')
      open(nu2,file=AF_out)

! ... read and write HEADER  and CLOSED

      read(nu1,'(a)') AS;    write(nu2,'(a)') AS
      read(nu1,'(a)') AS;    write(nu2,'(a)') AS

! ... read orbitals

      read(nu1,*) norb,ne,parity,in_corr,i_corr,j_corr

      Call alloc_orb_LS(norb)
      Allocate(iq_min(norb),iq_max(norb))

      read(nu1,'(50a4)') (ELF(i),i=1,norb)
      read(nu1,'(50i4)') (iq_min(i),i=1,norb)
      read(nu1,'(50i4)') (iq_max(i),i=1,norb)

      Do i=1,norb
       Call EL4_nlk(ELF(i),NEF(i),LEF(i),KEF(i))
      End do
      nwf=norb;  ncfg=0

      Close(nu1)

! ... generation of all possible configurations:

      Call Sum_conf(iq_min,iq_max,in_corr,i_corr,j_corr,nu2,mo)

      write(nu2,'(a)') '*'

      write(*,*)  ' ELECTRON.INP --> CONF.INP,   Nconf = ',ncfg


      End  ! end utility ZGENCONF


!----------------------------------------------------------------------
      Subroutine Sum_conf(iq_min,iq_max,in_corr,i_corr,j_corr,nu,mo)
!----------------------------------------------------------------------
!     exhaustion of possible configurations
!----------------------------------------------------------------------

      Use conf_LS,  only: ne
      Use orb_LS

      Integer :: iq_min(*),iq_max(*),in_corr,i_corr,j_corr,nu

      i1=1; i2=nwf; IEF=0 

      i=1;  IEF(i)=iq_max(i)

    1 ii=SUM(IEF(1:i))
      if(ii.eq.ne) then
        Call gen_conf(in_corr,i_corr,j_corr,nu,mo)
      elseif(ii.lt.ne.and.i.lt.i2) then
        i=i+1
        IEF(i)=iq_max(i)
        go to 1
      end if

    2 IEF(i)=IEF(i)-1
      if(IEF(i).lt.iq_min(i)) then
        if(i.eq.i1) Return
        i=i-1
        go to 2
      end if
      go to 1

      End  ! Subroutine Sum_conf


!====================================================================
      Subroutine Gen_conf(in_corr,i_corr,j_corr,nu,mo)
!====================================================================
!
!     generate and check 1 configuration with electron population
!     given by IEF in module configs
!
!--------------------------------------------------------------------

      USE conf_LS;  USE orb_LS

      Integer, intent(in) :: in_corr,i_corr,j_corr,nu

      no=0; n_corr=0
      Do i=1,nwf
       if(IEF(i).le.0) Cycle
       if(i.gt.in_corr) n_corr=n_corr+IEF(i)
       no=no+1
       nn(no)=NEF(i); ln(no)=LEF(i); kn(no)=KEF(i); iq(no)=IEF(i)
      End do

      if(n_corr.lt.i_corr) Return
      if(n_corr.gt.j_corr) Return
      if(no.gt.mo) Return

! ... check the total parity:

      k=0;  Do i=1,no; k=k+iq(i)*ln(i);  End do
      k=(-1)**k
      if(k.ne.parity) Return

! ... check the number of electrons in the shell

      Do i=1,no;  if(iq(i).gt.4*ln(i)+2) Return;  End do

! ... record configuration:

      Call Incode_c;   write(nu,'(a)') trim(CONFIG);   ncfg=ncfg+1

      End ! Subroutine Gen_conf


!======================================================================
      Subroutine create_electron_inp
!======================================================================
!     provide screen information about add_farm utility
!----------------------------------------------------------------------
       
      Implicit real(8) (a-h,o-z)

      nu=1;  open(nu,file='electron.inp') 

      write(nu,'(a)') &
'OI                                                                                ',&  
'  1s  2s  2p                                                                      ',&  
'   6   6   1   2   0   2   n_orbitals,  n_electrons, parity, n_ref, k_min, k_max  ',&  
'  3s  3p  4s  4p  4d  4f                                                          ',&  
'   0   0   0   0   0   0                                                          ',&  
'   2   6   2   2   2   2                                                          ',&  
'                                                                                  ',&  
'                                                                                  ',&
'n_orbitals  - number of orbitals in following list (up to 50 in a row)            ',&
'n_electrons - number of electrons above common core                               ',&
'parity      - -/+ 1                                                               ',&
'n_core      - reference orbitals for promotions                                   ',&
'k_min       - minumum promotion, e.g,   1 - single                                ',&
'k_max       - maximum promotion, e.g.,  2 - double                                '

      write(*,'(/a/)') 'example of electron.inp file was created  -  fill it out ! '

      Stop ' '

      End Subroutine create_electron_inp
                                        
