!=====================================================================
!     UTILITY  add_farm
!=====================================================================
!
!     accumulates data from different FARM calculations
!
!     Program first check the existance of the FARM output files
!
!       'farm.kma','farm.pha','farm.tma','farm.om'
!
!     and added the information if any to the FORMATTED files
!
!       'zarm.kma','zarm.pha','zarm.tma','zarm.om'
!
!     The FARM files then will be deleted
!
!     The program requires the 'target' file which contains all
!     scattering (target and channel) information
!     The 'target' file can be obtained from H.DAT file with
!     utility 'h_targb'
!
!     Call:   add_farm       (without arguments)
!
!---------------------------------------------------------------------

      Use target; Use channels

      Implicit real(8) (a-h,o-z)

      Character(80) :: AS

      Real(8), allocatable :: E(:)
      Real(8), allocatable :: OM(:,:)
      Integer, allocatable :: ichl(:),lchl(:)
      Real(8), allocatable :: matr(:),mati(:)

      Logical EX
      Character(20), Dimension(7) :: fname, zname

      Data fname/'farm.kma','farm.pha','farm.tma', &
                 'farm.sec','farm.om ','farm.top','farm.om'/
      Data zname/'zarm.kma','zarm.pha','zarm.tma', &
                 'zarm.sec','zarm.om ','zarm.top','zarm.om_par'/

      Call inf_add_farm

!----------------------------------------------------------------------
! ... target information:

      Inquire(file='target',EXIST=EX)
      if(.not.EX) then
       Stop  ' No target file: run H_targb first '
      else
       nut=1; Open(nut,file='target');
       Call R_target(nut)
       Call R_channels(nut)
       close(nut)
      end if

      E1 = etarg(1); ETARG = (ETARG-E1)*2
      ion = nz-nelc

!----------------------------------------------------------------------
! ... OMEGA   (itype=5):

      itype=5; inp=10+itype; iout=20+itype

      klsp=0; Call Read_iarg('klsp',klsp)
      if(klsp.gt.0) then
       itype=7; inp=10+itype; iout=20+itype
      end if

      Inquire (file=fname(itype),EXIST=EX)
      if(EX) then

      Open(inp,file=fname(itype),FORM='UNFORMATTED')
      Open(iout,file=zname(itype),POSITION='APPEND')

      read(inp) nt,ne
      if(Allocated(E)) Deallocate(E); Allocate(E(ne))
      read(inp) E(1:ne)
      if(Allocated(OM)) Deallocate(OM); Allocate(OM(nt,ne))
      read(inp) ((OM(it,ie),it=1,nt),ie=1,ne)

      Do ie=1,ne
       io=IOPEN(ntarg,E(ie),etarg)
       ns=io*(io-1)/2; if(ion.le.0) ns=io*(io+1)/2
       if(ns.eq.0) Cycle
       if(klsp.gt.0) &
       write(iout,'(F10.6,6i8)') E(ie),klsp,ns,io,io,ntarg,ntarg
       if(klsp.eq.0) &
       write(iout,'(F10.6,i8)') E(ie),ns
       write(iout,'(5D16.8)') (OM(it,ie),it=1,ns)
      End do

      Close(iout);  Close(inp,status='delete')
      Deallocate(E,OM)
      end if

!----------------------------------------------------------------------
! ... T-matrices     ( itype = 3)

      itype=3; inp=10+itype; iout=20+itype
      Inquire (file=fname(itype),EXIST=EX)
      if(EX) then

      Open(inp,file=fname(itype),FORM='UNFORMATTED')
      Open(iout,file=zname(itype),POSITION='APPEND')

      if(allocated(ichl)) Deallocate(ichl,lchl)
      Allocate(ichl(mch),lchl(mch))
      if(allocated(matr)) Deallocate(matr,mati)
      mdim=mch*(mch+1)/2; Allocate(matr(mdim),mati(mdim))

    3 read (inp,end=13) IS,IL,IP,ne,nchan,(ichl(i),lchl(i),i=1,nchan)
      IP = (-1)**IP

! ... find partial wave:

      ii=0
      Do i=1,nlsp
       if(IS.ne.ispar(i)) Cycle
       if(IL.ne. lpar(i)) Cycle
       if(IP.ne. ipar(i)) Cycle
       ii=i; Exit
      End do
      if(ii.eq.0) Stop ' Unknown symmetry in T-matrix'

      Do ie = 1,ne
       read(inp) ee,nopen,ntr,(matr(i),mati(i),i=1,ntr)
       write(iout,'(F10.6,3i8)') ee,nopen,ntr,ii
       write(iout,'(6D16.8)') (matr(i),mati(i),i=1,ntr)
      End do

      go to 3
   13 Continue

      Close(iout); Close(inp,status='delete')
      Deallocate(ichl,lchl,matr,mati)
      end if

!----------------------------------------------------------------------
! ... K-matrices     ( itype = 1)

      itype=1; inp=10+itype; iout=20+itype
      Inquire (file=fname(itype),EXIST=EX)
      if(EX) then

      if(allocated(ichl)) Deallocate(ichl,lchl)
      Allocate(ichl(mch),lchl(mch))
      if(allocated(matr)) Deallocate(matr)
      mdim=mch*(mch+1)/2; Allocate(matr(mdim))

      Open(inp,file=fname(itype),FORM='UNFORMATTED')
      Open(iout,file=zname(itype),POSITION='APPEND')


   1  read (inp,end=11) IS,IL,IP,ne,nchan,(ichl(i),lchl(i),i=1,nchan)
      IP = (-1)**IP

! ... find partial wave:

      ii=0
      Do i=1,nlsp
       if(IS.ne.ispar(i)) Cycle
       if(IL.ne. lpar(i)) Cycle
       if(IP.ne. ipar(i)) Cycle
       ii=i; Exit
      End do
      if(ii.eq.0) Stop ' Unknown symmetry in K_matrix'

      Do ie = 1,ne
       read(inp,end=11,err=11) ee,nopen,ntr,(matr(i),i=1,ntr)
       write(iout,'(F10.6,3i8)') ee,nopen,ntr,ii
       write(iout,'(5D16.8)') (matr(i),i=1,ntr)
      End do

      go to 1
   11 Continue

      Close(iout); Close(inp,status='delete')
      Deallocate(ichl,lchl,matr)
      end if

!----------------------------------------------------------------------
! ... phases         (itype = 2)

      itype=2; inp=10+itype; iout=20+itype
      Inquire (file=fname(itype),EXIST=EX)

      if(EX) then

      Open(inp,file=fname(itype))
      Open(iout,file=zname(itype),POSITION='APPEND')

   2  read(inp,'(a)',end=22) AS
      if(AS(21:21).ne.'.') go to 2
      read (AS,*) IS,IL,IP,ee,phase
      IP = (-1)**IP

! ... find partial wave:

      ii=0
      Do i=1,nlsp
       if(IS.ne. ispar(i)) Cycle
       if(IL.ne. lpar(i)) Cycle
       if(IP.ne. ipar(i)) Cycle
       ii=i; Exit
      End do
      if(ii.eq.0) Stop ' Unknown symmetry in K_matrix'

      write(iout,'(2F16.8,4i5)') ee,phase,ii,IS,IL,IP

      go to 2
      end if
   22 Close(iout); Close(inp,status='delete')

!----------------------------------------------------------------------

      End  ! UTILITY  add_farm



!======================================================================
      Integer Function Iopen(n,E,ETARG)
!======================================================================
!     number of open target states for given energy E
!----------------------------------------------------------------------

      Implicit none

      Integer,intent(in) :: n
      Real(8),intent(in) :: E
      Real(8),intent(in), Dimension (n) :: Etarg
      Integer :: i

      Iopen=0
      Do i=1,n
       if(E.gt.ETARG(i)) Iopen=i
      End do

      End Function Iopen


!======================================================================
      Subroutine inf_add_farm
!======================================================================
!     provide screen information about add_farm utility
!----------------------------------------------------------------------

      Character :: A

      Call get_command_argument(1,A)

      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                ',&
'     ADD_FARM accumulates data from different FARM calculations ',&
'                                                                ',&
'     Program first check the existance of the FARM output files ',&
'                                                                ',&
'       farm.kma, farm.pha, farm.tma, farm.om                    ',&
'                                                                ',&
'     and added the information if any to the FORMATTED files    ',&
'                                                                ',&
'       zarm.kma, zarm.pha, zarm.tma, zarm.om                    ',&
'                                                                ',&
'     The FARM files then will be deleted                        ',&
'                                                                ',&
'     The program requires the target file which contains all    ',&
'     scattering (target and channel) information                ',&
'     The target file can be obtained from H.DAT file with       ',&
'     utility h_targb                                            ',&
'                                                                ',&
'     Call as:   add_farm    [klsp=..]                           ',&
'                                                                ',&
'     klsp - index of partial wave, if calculations for one partial wave only,',&
'            in this case output is recorded in zarm.om_par file ',&
'                                                                ',&
'                                                                '
      Stop ' '

      End Subroutine inf_add_farm

