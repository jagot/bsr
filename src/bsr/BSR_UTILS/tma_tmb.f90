!======================================================================
!     UTILITY  tma_tmb
!======================================================================
!
!     zarm.tma  --> zarm.tmb
!
!     provide new (reduced) format for t-matrix file
!
!     Call:  tma_tmb  [tma=.. tmb=.. ]
!
!     np -  number of physical state which will be kept in t-matrix file
!     ni -  number of states with all excitations possible
!======================================================================

      Use target; Use channels

      Implicit real(8) (a-h,o-z)
      Real(8), Allocatable ::  tmar(:),tmai(:)

! ... files:

      Character(80) :: targ  = 'target';     Integer :: nut = 1
      Character(80) :: tma   = 'zarm.tma';   Integer :: nu1 = 2
      Character(80) :: tmb   = 'zarm.tmb';   Integer :: nu2 = 3

      Integer :: np = 0,  ni = 0

      Character(80) :: A

      Call get_command_argument(1,A)
      if(A.eq.'?'.or. A.eq.'!') then

      write(*,'(a)') &
'                                                                          ',&
'     zarm.tma  --> zarm.tmb                                               ',&
'                                                                          ',&
'     provides new (reduced) format for t-matrix file                      ',&
'                                                                          ',&
'     Call as:  tma_tmb  [tma=zarm.tma  tmb=zarm.tmb  np=ntarg  ni=ntarg]  ',&
'                                                                          ',&
'     np,ni -  if differ from ntarg, should be given in "target" file      ',&
'                                                                          ',&
'     default values are indicated in the example                          ',&
'                                                                          '

      STOP ' '
      end if

!----------------------------------------------------------------------
! ... target information:

      Call Check_file(targ)
      Open(nut,file=targ)
      Call R_target(nut)
      Call R_channels(nut)

      np=ntarg; Call Read_ipar(nut,'np',np)
      ni=ntarg; Call Read_ipar(nut,'ni',ni)

      Close(nut)

      mdim=mch*(mch+1)/2
      Allocate(tmar(mdim),tmai(mdim))

! ... re-define input-output files if needed:

      Call Read_aarg('tma',tma)
      Call Read_aarg('tmb',tmb)

      write(*,'(/a,a)')  'tma =',tma
      write(*,'( a,a)')  'tmb =',tmb

      write(*,'(/a)') 'read the T-matrix from tma-file and added to tmb-file (in new style)'

!----------------------------------------------------------------------
! ... parameters:

      Call Read_iarg('np',np);   if(np.le.0.or.np.gt.ntarg) np = ntarg
      Call Read_iarg('ni',ni);   if(ni.gt.ntarg) ni = ntarg

      write(*,'(/a,i5,a)')  'np =',np,' - number of physical states'
      write(*,'( a,i5,a)')  'ni =',ni,' - number of ionized states'

!----------------------------------------------------------------------

      Call Check_file(tma)
      open(nu1,file=tma)
      open(nu2,file=tmb,position='APPEND')

    1 read(nu1,*,end=2) ee,nopen,ntr,ilsp

      if(ntr.gt.mdim) Stop 'ntr > mdim in T-matrix file'
      if(ilsp.gt.nlsp) Stop 'ilsp > nlsp in T-matrix file'
      if(nopen.gt.nch(ilsp)) Stop 'nopen > nch(ilsp)'

      read(nu1,*) (tmar(i),tmai(i),i=1,ntr)

      i = Jopen(ee,ilsp)
      if(nopen.gt.i) then
       write(81,*) 'i,nopen,ilsp,e',i,nopen,ilsp,ee
       go to 1
      end if

      kp = 0;  nj = 0
      Do ich = 1,nch(ilsp)
       if(iptar(ilsp,ich).le.np) kp=ich
       if(iptar(ilsp,ich).le.ni) nj=ich
      End do

      if(kp.eq.0) go to 1

      if(kp.gt.nopen) kp=nopen

      write(nu2,'(F10.6,6i6,a)') ee,nopen,kp,ilsp,np,ni,nj ,'   ee,nopen,kp,ilsp,np,ni,nj'

      write(nu2,'(6D16.8)') &
         ((tmar(i*(i-1)/2+j),tmai(i*(i-1)/2+j),j=1,i),i=1,kp)

      if(nopen.gt.kp.and.nj.gt.0)   write(nu2,'(6D16.8)') &
         ((tmar(i*(i-1)/2+j),tmai(i*(i-1)/2+j),j=1,nj),i=kp+1,nopen)

      go to 1
    2 Continue

      Close(nu1)
      Close(nu2)

      End  !   tma_tmb










