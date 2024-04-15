!=====================================================================
!     UTILITY  convert_target
!=====================================================================
!     target, threshoulds -> target.exp
!
!     Call as:  indicated are the optional (default) paremeters:
!
!     convert_target  [ targ1=target.theory
!                       targ2=target.exp
!                       thresholds=thereholds ]
!---------------------------------------------------------------------

      Use target; Use channels

      Implicit real(8) (A-H,O-Z)

      Real(8), allocatable :: E_exp(:)
      Integer, allocatable :: ip_exp(:), jp_exp(:), npch(:), nptar(:)
      Character(72) :: title, AF

      Integer :: nu1=1;  Character(80) :: AF1  = 'target.theory'
      Integer :: nu2=2;  Character(80) :: AF2  = 'target.exp'
      Integer :: nu3=3;  Character(80) :: AF3  = 'thresholds'

!----------------------------------------------------------------------
      iarg = COMMAND_ARGUMENT_COUNT()
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)

      if(AF.eq.'?') then
        write(*,'(/a)') 'change target file to fit new experimental state energies'
        write(*,'(/a)') 'given in file thresholds'
        write(*,'(/a)') 'default:   target.theory + thrsholds => target.exp'
        write(*,'(/a)') 'change default names as:'
        write(*,'(/a)') 'convert_target  [targ1=...  targ2=... thresholds=...]'
        Stop
      end if



      Call Read_aarg('targ1',AF1)
      Call Read_aarg('targ2',AF2)
      Call Read_aarg('thresholds',AF3)

!----------------------------------------------------------------------
! ... original target and channel information:

      Open(nu1,file=AF1,status='OLD')
      Call R_target (nu1)
      Call R_channels(nu1)
      rewind(nu1); read(nu1,'(a)') title
      Close(nu1)

      Z = nz
      Call Conv_au (Z,0.d0,au_cm,au_eV,0)

!----------------------------------------------------------------------
! ... new thresholds:

      Open(nu3,file=AF3,status='OLD')

      Allocate(E_exp(ntarg))
      Do i=1,ntarg; read(nu3,*) E_exp(i); End do
      Allocate(ip_exp(ntarg),jp_exp(ntarg))
      Call SORTR(ntarg,E_exp,ip_exp)
      iiexp=0
      Do i=1,ntarg; if(ip_exp(i).ne.i) iiexp=1; End do
      etarg = E_exp

      if(iiexp.eq.0) write(*,*) 'target states order is not changed'
      if(iiexp.eq.1) write(*,*) 'target states order is changed !!!'

!----------------------------------------------------------------------
! ... new target:

      Open(nu2,file=AF2);  nut=nu2

      write(nut,'(a)') TITLE
      write(nut,'(72(''-''))')
      write(nut,'(a,a2,4x,a)') &
                'coupling = ',coupling, '!   coupling scheme'
      write(nut,'(a,i5,T18,a)') &
                'nz     =',nz,          '!   nuclear charge'
      write(nut,'(a,i5,T18,a)') &
                'nelc   =',nelc,        '!   number of electrons'
      write(nut,'(72(''-''))')
      write(nut,'(a,i5,T18,a)') &
                'ntarg  =',ntarg,       '!   number of target states'
      write(nut,'(72(''-''))')
      Do j=1,ntarg; i=ip_exp(j); jp_exp(i) = j
       E_Ry = (etarg(i)-etarg(1))*2
       E_eV = (etarg(i)-etarg(1))*au_eV
       write(nut,'(a20,1x,a10,1x,3i4,f18.8,f10.6,f10.3)') BFT(i),AFT(i), &
        ltarg(i),istarg(i),iptarg(i),etarg(i), E_Ry, E_eV
      End do
      write(nut,'(72(''-''))')
      write(nut,'(a,i5,T18,a)') 'nct    =',nct, &
       '!   total number of target configurations'
      write(nut,'(a,i5,T18,a)') 'nwt    =',nwt, &
       '!   total number of target orbitals'

      write(nut,'(72(''-''))')
      write(nut,'(a,i5,T18,a)') &
           'nlsp   =',nlsp, '!   number of partial waves'
      write(nut,'(72(''-''))')

      AFP = 'no';  BFP='no'
      Do i = 1,nlsp
       if(AFP(i).ne.'no'.or.BFP(i).ne.'no') then
         write(nut,'(i3.3,3i5,3x,a20,1x,a10,2i5)') &
         i,lpar(i),ispar(i),ipar(i),AFP(i),BFP(i),ncp(i),nwp(i)
       else
         write(nut,'(i3.3,3i5)') i,lpar(i),ispar(i),ipar(i)
       end if
      End do
      write(nut,'(72(''-''))')

      write(nut,'(a,i5)') 'channels:'
      write(nut,'(80(''-''))')

      Do ilsp = 1,nlsp
       write(nut,'(i3,a,i5.3,a,i6,a,2i10)') ilsp,'.',ilsp,' nch = ',nch(ilsp), &
        ' nc = ',ipconf(ilsp,nch(ilsp)),ncp(ilsp)

       allocate( npch(nch(ilsp)), nptar(nch(ilsp)) )

       ! ... new channel order:

       k = 0
       Do i=1,ntarg; it=ip_exp(i)
        Do ich=1,nch(ilsp); if(iptar(ilsp,ich).ne.it) Cycle
         k=k+1; npch(k) = ich; nptar(k) = i
        End do
       End do
       if(k.ne.nch(ilsp)) Stop 'Problems with new channel order'

       Do j = 1,nch(ilsp); i = npch(j)
        write(nut,'(a4,2x,3i6,2i8)') &
         ELC(ilsp,i),lch(ilsp,i),nptar(j),i,ipconf(ilsp,i),jkch(ilsp,i)
       End do

       write(nut,'(80(''-''))')

       deallocate(npch,nptar)

      End do


      End ! program
