!=====================================================================
!     UTILITY  add_stgfb
!=====================================================================
!     accumulates results STGF calculations (v.4.7.1 and pstgf format)
!
!     OMEGA, KMAT.DAT + target  --> zarm.kma, zarm.om, zarm.om_par
!
!     Call:   add_stgf   [klsp=..]
!
!---------------------------------------------------------------------
      Use target; Use channels

      Implicit real(8) (a-h,o-z)

      Character(80) :: AS
      Real(8), allocatable :: kmat(:,:), tmat(:,:), om(:,:), omb(:)

      Integer, parameter :: nproc = 9999

! ... files:

      Integer :: inp=1;   Character(80) :: tname   = 'TMAT.DAT'
                          Character(80) :: kname   = 'KMAT.DAT'
                          Character(80) :: oname   = 'OMEGA'
      Integer :: out=2;   Character(80) :: ztname  = 'zarm.tma'
                          Character(80) :: btname  = 'zarm.tmb'
                          Character(80) :: zkname  = 'zarm.kma'
                          Character(80) :: zoname  = 'zarm.om'
                          Character(80) :: boname  = 'zarm.omb'
                          Character(80) :: zpname  = 'zarm.om_par'
                          Character(80) :: bpname  = 'zarm.omb_par'
      Integer :: nut=3;   Character(80) :: targ    = 'target'
      Integer :: pri=6;   Character(80) :: AF_log  = 'add_stgf.log'
                          Character(80) :: fname

      Call inf_add_stgf
      open(pri,file=AF_log)

      Call CPU_time(t1)

!----------------------------------------------------------------------
! ... target information:

      Call Check_file(targ)
      Open(nut,file=targ)
      Call R_target(nut)
      Call R_channels(nut)

      np=ntarg; Call Read_ipar(nut,'np',np)
      ni=ntarg; Call Read_ipar(nut,'ni',ni)

      Call Read_iarg('np',np)
      Call Read_iarg('ni',ni)

      Close(nut)

      E1 = etarg(1); ETARG = (ETARG-E1)*2
      ion = nz-nelc;  zion=ion; if(ion.eq.0) zion=1.d0

      Allocate(kmat(mch,mch),tmat(mch,mch),OM(ntarg,ntarg), &
               OMB(ntarg*(ntarg+1)/2))

!----------------------------------------------------------------------
! ... read K-matrix:

      Open(out,file=zkname,POSITION='APPEND')

      Do ii = -1,nproc

       if(ii.eq.-1) fname = kname
       if(ii.gt.-1) write(fname,'(a,i4.4)') trim(kname),ii
       i = Icheck_file(fname)
       if(i.eq.0) Cycle

       Open(inp,file=fname)

    1 Continue
      if(coupling.eq.'LS') then
       read(inp,'(3i5,e14.5)',end=2) nopen,IST,ILT,EK
       IPT = 1; if(IST.lt.0) IPT=-1
       IST = iabs(IST)
      else
       read(inp,'(3i5,e14.5)',end=2) nopen,ILT,IPT,EK
       IPT = IPT + 1;  if(IPT.eq.2) IPT=-1
       IST = 0
      end if

      if(nopen.gt.mch) write(pri,*) 'Stop: open > mch'
      if(nopen.gt.mch) Stop 'nopen > mch'
      Do i=1,nopen; read(inp,*) lt,lc; End do
      if(nopen.le.0.and.ii.eq.-1) read(inp,*) lt
      if(nopen.le.0) go to 1

      ilsp = 0
      Do i=1,nlsp
       if(IST.ne.ispar(i)) Cycle
       if(ILT.ne.lpar (i)) Cycle
       if(IPT.ne.ipar (i)) Cycle
       ilsp =i; Exit
      End do

      if(ilsp.eq.0) then
       write(*,*) 'IST,ILT,IPT=',IST,ILT,IPT
       write(pri,*) 'IST,ILT,IPT=',IST,ILT,IPT
       write(pri,*) 'Stop: unknown symmetry in KMAT.DAT'
       Stop 'unknown symmetry in KMAT.DAT'
      end if

      Do i=1,nopen; Do j=1,nopen
       read(inp,*) kmat(i,j)
      End do; End do

! ... save K-matrix:

      ntr = nopen*(nopen+1)/2

      write(out,'(F10.6,3i10,a)') ek,nopen,ntr,ilsp, '  ek,nopen,ntr,ilsp'
      write(out,'(5d15.6)') ((KMAT(i,j),i=1,j),j=1,nopen)
      write(pri,'(a,5x,a,F15.8,a,i5,a,i10,a,i4)') trim(fname), &
       'ek =',ek,'  nopen=',nopen,'  ntr=',ntr,'  ilsp=',ilsp

      go to 1
    2 Continue
      Close(inp,status='delete')

      End do   ! over ii -> file index

!----------------------------------------------------------------------
! ... read T-matrix:

      if(np.eq.ntarg) Open(out,file=ztname,POSITION='APPEND')
      if(np.lt.ntarg) Open(out,file=btname,POSITION='APPEND')

      Do ii = -1,nproc
       if(ii.eq.-1) fname = tname
       if(ii.gt.-1) write(fname,'(a,i4.4)') trim(tname),ii
       i = Icheck_file(fname)
       if(ii.eq.-1.and.i.eq.0) Cycle
       if(ii.gt.-1.and.i.eq.0) Cycle  ! Exit

       Open(inp,file=fname)

   10 Continue
      if(coupling.eq.'LS') then
       read(inp,'(3i5,e14.5)',end=20) nopen,IST,ILT,EK
       IPT = 1; if(IST.lt.0) IPT=-1;  IST = iabs(IST)
      else
       read(inp,'(3i5,e14.5)',end=20) nopen,ILT,IPT,EK
       IPT = IPT + 1;  if(IPT.eq.2) IPT=-1;  IST = 0
      end if

      if(nopen.gt.mch) write(pri,*) 'Stop: open > mch'
      if(nopen.gt.mch) Stop 'nopen > mch'
      Do i=1,nopen; read(inp,*) lt; End do
      if(nopen.le.0) then; read(inp,*) lt;  go to 10; end if

      ilsp = 0
      Do i=1,nlsp
       if(IST.ne.ispar(i)) Cycle
       if(ILT.ne.lpar (i)) Cycle
       if(IPT.ne.ipar (i)) Cycle
       ilsp =i; Exit
      End do

      if(ilsp.eq.0) then
       write(*,*) 'IST,ILT,IPT=',IST,ILT,IPT
       write(pri,*) 'IST,ILT,IPT=',IST,ILT,IPT
       write(pri,*) 'Stop: unknown symmetry in KMAT.DAT'
       Stop 'unknown symmetry in KMAT.DAT'
      end if

      Do i=1,nopen; Do j=1,nopen
       read(inp,*) kmat(i,j),tmat(i,j)
      End do; End do

! ... save T-matrix:

      ntr = nopen*(nopen+1)/2

      if(np.eq.ntarg) then

       write(out,'(F10.6,3i10,a)') ek,nopen,ntr,ilsp, '  ek,nopen,ntr,ilsp'
       write(out,'(6d15.6)') ((KMAT(i,j),TMAT(i,j),i=1,j),j=1,nopen)
       write(pri,'(a,5x,a,F15.8,a,i5,a,i10,a,i4)') trim(fname), &
        'ek =',ek,'  nopen=',nopen,'  ntr=',ntr,'  ilsp=',ilsp

      else

       i = Jopen(ek,ilsp)
       if(nopen.gt.i) then
        write(pri,*) 'WARNING: i < nopen, ilsp, e',i,nopen,ilsp,ee
        go to 10
       end if

       kp = 0;  nj = 0
       Do ich = 1,nch(ilsp)
        if(iptar(ilsp,ich).le.np) kp=ich
        if(iptar(ilsp,ich).le.ni) nj=ich
       End do

       if(kp.eq.0) go to 10

       if(kp.gt.nopen) kp=nopen

       write(out,'(F10.6,6i6,a)') ek,nopen,kp,ilsp,np,ni,nj , &
                              '   ee,nopen,kp,ilsp,np,ni,nj'

       write(out,'(6D16.8)')  ((KMAT(i,j),TMAT(i,j),i=1,j),j=1,kp)  ! lower part

       if(nopen.gt.kp.and.nj.gt.0)   write(out,'(6D16.8)') &
          ((KMAT(i,j),TMAT(i,j),j=1,nj),i=kp+1,nopen)

      end if

      go to 10
   20 Continue
      Close(inp,status='delete')

      End do   ! over ii -> file index

!----------------------------------------------------------------------
! ... read OMEGA:

      klsp=0; Call Read_iarg('klsp',klsp)

      if(np.eq.ntarg.and.klsp.eq.0)  fname=zoname
      if(np.eq.ntarg.and.klsp.gt.0)  fname=zpname

      if(np.ne.ntarg.and.klsp.eq.0)  fname=boname
      if(np.ne.ntarg.and.klsp.gt.0)  fname=bpname

      Open(out,file=fname,POSITION='APPEND')

      Do ii = -1,nproc

       write(fname,'(a,i4.4)') trim(oname),ii
       if(ii.eq.-1) fname = trim(oname)
       i = Icheck_file(fname)
       if(i.eq.0) Cycle

!       if(i.eq.0.and.ii.eq.-1) Cycle
!       if(i.eq.0.and.ii.gt.-1) Exit

       Open(inp,file=fname)

       READ(inp,*) NZED,NELC1
       if(nz.ne.NZED) Stop ' nz <> NZED'
       if(nelc.ne.NELC1) Stop ' nelc ?'

       READ(inp,*) NAST,ne,NOMWRT
       if(ntarg.ne.NAST) Stop ' ntrag <> NAST'
       if(NOMWRT.eq.0) Stop 'NOMWRT = 0 '

       READ(inp,*) (ISAT,LAT,I=1,NAST)
       READ(inp,*) (ENAT,I=1,NAST)

       Do ie = 1,ne

       om = 0.d0
       if(nz.eq.nelc) then
        read(inp,*) ek,((OM(i,j),i=j,NAST),j=1,NAST)
       else
        read(inp,*) ek,((OM(i,j),i=j+1,NAST),j=1,NAST-1)
       end if

       Do j=1,NAST; Do i=j,NAST
        OM(j,i) = OM(i,j)
       End do; End do

       ek = ek*zion*zion

       io=IOPEN(ntarg,ek,etarg)

       if(np.eq.ntarg) then

        if(nz.eq.nelc) then
         ns=io*(io+1)/2
         if(klsp.gt.0) write(out,'(F10.6,6i8)') ek,klsp,ns,io,io,ntarg,ntarg
         if(klsp.eq.0) write(out,'(F10.6, i8)') ek,ns
         write(out,'(5D16.8)') ((OM(i,j),j=1,i),i=1,io)
        else
         ns=io*(io-1)/2
         if(klsp.gt.0) write(out,'(F10.6,6i8)') ek,klsp,ns,io,io,ntarg,ntarg
         if(klsp.eq.0) write(out,'(F10.6, i8)') ek,ns
         write(out,'(5D16.8)') ((OM(i,j),j=1,i-1),i=2,io)
        end if

       else

        jop = io

        iop = jop; if(iop.gt.np) iop=np
        if(ion.eq.0) ntr = iop*(iop+1)/2
        if(ion.ne.0) ntr = iop*(iop-1)/2
        ktr = ntr; if(jop.gt.iop) ktr = ntr + ni*(jop-iop)

        omb = 0.d0
        Do j = 1,iop; Do i=1,j
         in = Index_TR(ion,i,j,np,ni)
         if(in.eq.0) Cycle
         omb(in)=OM(i,j)
        End do; End do

        Do j=iop+1,jop; Do i=1,ni
         in = Index_TR(ion,i,j,np,ni)
         if(in.eq.0) Cycle
         omb(in)=OM(i,j)
        End do; End do

        if(klsp.eq.0) &
         write(out,'(F10.6,5i8,a)') ek,ktr,iop,jop,np,ni,'   e,ntr,iopen,jopen,np,ni'
        if(klsp.ne.0) &
         write(out,'(F10.6,6i8,a)') ek,klsp,ktr,iop,jop,np,ni, &
                                 '  e,ilsp,nom,iopen,jopen,np,ni'
        write(out,'(5D16.8)') (omb(i),i=1,ktr)

        write(pri,'(a,5x,a,F15.8,a,i4,a,i10,a,i10)') trim(fname), &
          'ek =',ek,'  nopen=',io,'  ntr=',ns,'  ilsp=',klsp

       end if

       End do  ! ie

       Close(inp,status='delete')

       End do  ! over ii

      Call CPU_time(t2)
      write(pri,'(/a,f10.1,a)') 'time: ',(t2-t1)/60,' min'

      End  ! UTILITY  add_stgf


!======================================================================
      Subroutine inf_add_stgf
!======================================================================
!     provide screen information about add_stgf utility
!----------------------------------------------------------------------

      Character :: A = ' '

      Call get_command_argument(1,A)

      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                           ',&
'     add_stgf accumulates results after STGF (PSTGF) runs                  ',&
'                                                                           ',&
'     OMEGA, KMAT.DAT, TMAT.DAT + target  --> zarm.om, zarm.kma, zarm.tma   ',&
'     (if np /= ntarg)                    --> zarm.omb, zarm.kma, zarm.tmb ',&
'                                                                           ',&
'     Call as:   add_stgfb   [klsp=..]                                      ',&
'                                                                           ',&
'     klsp - index of partial wave, if calculations for one                 ',&
'     partial wave only, in this case output is recorded in                 ',&
'     zarm.om_par|zarm.omb_par file                                         ',&
'                                                                           ',&
'     np   - number of physical states                                      ',&
'                                                                           '

      Stop ' '

      End Subroutine inf_add_stgf

