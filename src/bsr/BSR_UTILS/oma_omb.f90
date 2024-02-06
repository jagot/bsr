!=====================================================================
!     UTILITY  oma_omb
!=====================================================================
!
!     Change fortmat from zarm.omb to original zarm.om
!
!     Arguments: omb  -   zarm.omb file
!                om   -   zarm.om  file
!
!---------------------------------------------------------------------
      Use target

      Implicit real(8) (A-H,O-Z)

      Real(8), Allocatable :: om(:),y(:)

      Character(20)  :: AT

      Real(8) :: Ry = 13.6057

      Integer :: nut=1;  Character(80) :: AF_t   = 'target'
      Integer :: nub=2;  Character(80) :: AF_omb = 'zarm.omb'
      Integer :: nua=3;  Character(80) :: AF_oma = 'zarm.om'

      Call inf_sub

!----------------------------------------------------------------------
! ... target information:

      Call Check_file(AF_t)
      Open(nut,file=AF_t)
      Call R_target(nut)
      np=ntarg; Call Read_ipar(nut,'np',np)
      ni=ntarg; Call Read_ipar(nut,'ni',ni)
      close(nut)

      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion

      E1 = etarg(1); etarg = (etarg-E1)*2.0

      Z = nz;  AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,0)
      Ry = au_eV/2.d0

! ... define files:

      Call Read_aarg('oma',AF_oma)
      Call Check_file(AF_oma)
      open(nua,file=AF_oma)

      Call Read_aarg('omb',AF_omb)
      Call Check_file(AF_omb)
      open(nub,file=AF_omb)

      mdim = ntarg*(ntarg+1)/2
      Allocate(y(mdim),om(mdim))

!----------------------------------------------------------------------
! ... transfer for each energies:

    1 read(nua,*,end=2) E,ns, (y(i),i=1,ns)

      om = 0.d0

      jop = Iopen(ntarg,E,ETARG)

      if(ion.eq.0.and.jop*(jop+1)/2.ne.ns.or. &
       ion.ne.0.and.jop*(jop-1)/2.ne.ns) then
        write(*,*) E,ns, jop
        go to 1
      end if

      iop = jop; if(iop.gt.np) iop=np
      if(ion.eq.0) ntr = iop*(iop+1)/2
      if(ion.ne.0) ntr = iop*(iop-1)/2
      ktr = ntr; if(jop.gt.iop) ktr = ntr + ni*(jop-iop)

      om(1:ntr)=y(1:ntr)

      if(jop.gt.iop) then
       if(iop.ne.np)  Stop 'iop <> np ?'
       Do i=iop+1,jop
        Do j=1,ni
         itr = i*(i-1)/2 + j
         if(itr.gt.ns) Stop 'itr > ns'
         jtr = ntr + (i-iop-1)*ni + j
         if(jtr.gt.ktr) Stop 'jtr > ktr'
         om(jtr)=y(itr)
        End do
       End do
      end if

      write(nub,'(F10.6,5i8,a)')  e,ktr,iop,jop,np,ni, &
                              '   e,ntr,iopen,jopen,np,ni'
      write(nub,'(5D16.8)') (om(i),i=1,ktr)

      go to 1
    2 Continue

      End  ! UTILITY  omb_oma


!======================================================================
      Subroutine inf_sub
!======================================================================
!     provide screen information about sec_om utility
!----------------------------------------------------------------------

      Character :: A

      iarg = command_argument_count()
      if(iarg.eq.0) Return
      Call GET_COMMAND_ARGUMENT(1,A)
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                ',&
'     utility omb_oma changes format from zarm.om to zarm.omb    ',&
'                                                                ',&
'          zarm.om    =>   zarm.omb                              ',&
'                                                                ',&
'     Call as:  oma_omb  [ oma=..  omb=..]                       ',&
'                                                                ',&
'     Parameters oma and omb may be used to redefine the input/output files ',&
'                                                                '
      Stop ' '

      End Subroutine inf_sub

