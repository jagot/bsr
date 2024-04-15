!======================================================================
!     UTILITY      S E C _ T O P _TM
!
!     zarm.tma (zarm.tmb) --> zarm.oma_top (zarm.tmb_top)
!
!======================================================================
!     genearte top-up omegas for all included transitions
!     (based on the information in 'zarm.tmb)
!
!     Call as:  sec_top  [tmb= top= ek1= ek2= tail= x= method= ...]
!---------------------------------------------------------------------

      Use target; Use channels

      Implicit real(8) (a-h,o-z)

      Real(8), allocatable ::  fl(:), om(:), e(:), fail(:),coefa(:),coefb(:)
      Real(8), allocatable ::  om_sum(:,:), om_top(:,:), fll(:,:)
      Real(8), allocatable ::  fom(:,:,:), tmatr(:), tmati(:), OS(:,:)

      Integer, allocatable ::  iop(:,:,:), met(:), ic(:), jc(:)

      Integer :: ke = 20000  !  initial number of energies

      Real(8) :: eps_tail = 0.001
      Real(8) :: eps_x    = 0.025

      Integer :: np=0, ni=0, method=1, mem=0

      Character :: form
      Character(80) ::  AS

      Real(8), external :: Fdip0, IZ0_lamda
      Real(8) :: kap1, kap2
!----------------------------------------------------------------------
! ... files:

      Character(80) :: AF, label=' '
      Integer :: nup=11;  Character(80) :: targ  = 'target'
      Integer :: nut=12;  Character(80) :: tm
                          Character(80) :: tma   = 'zarm.tma'
                          Character(80) :: tmb   = 'zarm.tmb'
      Integer :: nuo=14;  Character(80) :: top   = 'zarm.omb_top'
      Integer :: nuq=15;  Character(80) :: oms   = 'zarm.omb'
      Integer :: pri= 6;  Character(80) :: out   = 'sec_top_omb.log'
      Integer :: nuc=17;  Character(80) :: ccc   = 'sec_top_coef_fail'
      Integer :: nub=18;  Character(80) :: bad   = 'bad_energies'
      Integer :: nus=19;  Character(80) :: AF_OS = 's_values'

      Integer :: nua=99;  ! scratch file

      Call CPU_time(t1)

!----------------------------------------------------------------------

      Call Inf_top_omb

      Open(pri,file=out)

!----------------------------------------------------------------------
! ... target information:

      Call Check_file(targ)
      Open(nup,file=targ)
      Call R_target(nup)
      Call R_channels(nup)
      np=ntarg; Call Read_ipar(nup,'np',np)
      ni=ntarg; Call Read_ipar(nup,'ni',ni)
      Close(nup)
      e1=etarg(1); Do i=1,ntarg; etarg(i)=(etarg(i)-e1)*2; End do
      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion

      Z = nz
      AWT = 0.d0;  Call Read_rarg('AWT',AWT)
      Call Conv_au (Z,AWT,au_cm,au_eV,0)
      Ry = au_eV/2.d0

!-----------------------------------------------------------------------
! ... define other parameters:

      ek1 = 0.d0; Call Read_rarg('ek1',ek1)
      ek2 = 0.d0; Call Read_rarg('ek2',ek2)
      ekk = 0.d0; Call Read_rarg('ekk',ekk)

      eps_tail=0.001; Call Read_rarg('tail',eps_tail)
      eps_x=0.025;    Call Read_rarg('x',eps_x)
      method = 1;     Call Read_iarg('method',method)

      jtr1=0; Call Read_iarg('jtr1',jtr1)
      jtr2=0; Call Read_iarg('jtr2',jtr2)

      Call Read_aarg('label',label)

      write(*,'(/a,f10.3,a)') 'x     =',eps_x,' - minimum geometric series factor'
      write(*,'( a,f10.3,a)') 'tail  =',eps_tail,' - correction to be worried'

      if(ek1.ne.0.d0) &
      write(*,'(/a,f10.6,a)') 'ek1   =',ek1,' - minimum electron energy allowed'
      if(ek2.ne.0.d0) &
      write(*,'( a,f10.6,a)') 'ek1   =',ek2,' - maximum electron energy allowed'
      if(ekk.ne.0.d0) &
      write(*,'( a,f10.6,a)') 'ekk   =',ekk,' - min. electron energy for extrapolation'

      if(jtr1.ne.0.d0) &
      write(*,'(/a,i10,a)') 'jtr1    =',jtr1,' - initial state'
      if(jtr2.ne.0.d0) &
      write(*,'( a,i10,a)') 'jtr2    =',jtr2,' - final state'

!-----------------------------------------------------------------------
! ... get s-values:

      open(nus,file=AF_OS)
      read(nus,*) nsvalues
      Allocate(OS(ntarg,ntarg)); OS =0.d0; mem = mem + 2*ntarg*ntarg

      Do n=1,nsvalues
       read(nus,*) i,j,OS(i,j); OS(j,i) = OS(i,j)
      End do

!----------------------------------------------------------------------
! ... find energies:

      form = 'b'
      if(Icheck_file(tmb).gt.0) then
        tm=tmb; form='b'
      elseif(Icheck_file(tma).gt.0) then
        tm=tma; form='a'
      end if

      Call Read_aarg('tm',tm)
      Call Read_aarg('form',form)
      Call Check_file(tm)
      Open(nut,file=tm)

      me = ke; Allocate(e(me)); mem=mem+2*me

      ne=0; mom=0
    1 if(form.eq.'a') then
       read(nut,*,end=2) e1,nopen,nom,ilsp; i1=np; i2=ni
       read(nut,*) (S,S,i=1,nom)
      else
       read(nut,*,end=2) e1,nopen,kopen,ilsp,i1,i2,nj
       ntr = kopen*(kopen+1)/2
       read(nut,*) (S,S,i=1,ntr)
       ktr = (nopen-kopen)*nj
       if(ktr.gt.0) read(nut,*) (S,S,i=ntr+1,ntr+ktr)
       nom = ntr+ktr
      end if

      if(ek1.gt.0.d0.and.e1.lt.ek1) go to 1
      if(ek2.gt.0.d0.and.e1.gt.ek2) go to 1

      if(np.ne.i1) Stop 'different np'
      if(ni.ne.i2) Stop 'different ni'

      if(nom.gt.mom) mom=nom

      ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

      if(ie.eq.0) then; ne=ne+1; e(ne)=e1;  end if

      if(ne.eq.me) then
       open(nua,form='UNFORMATTED',status='SCRATCH')
       rewind(nua);   write(nua) (e(i),i=1,ne)
       Deallocate(e); me=me+ke; Allocate(e(me))
       rewind(nua);   read(nua) (e(i),i=1,ne)
      end if

      go to 1
    2 write(*,'(/a,i5,a)')    'ne    = ',ne,' - number of energies'
      if(ne.eq.0) Stop 'nothing to do !'

      Call Rsort(ne,e)

      write(*,'(/a,f15.6,a)') 'e(1)  = ', e(1),  ' - min. electron energy'
      write(*,'( a,f15.6,a)') 'e(ne) = ', e(ne), ' - min. electron energy'

      write(*,'(/a,i5,a)')    'mom   = ', mom,   ' - maximum matrix dimension'

      write(*,'(/a,i5,a)')    'np    = ', np,    ' - number of physical states'
      write(*,'( a,i5,a)')    'ni    = ', ni,    ' - number of "ionization" states'

!----------------------------------------------------------------------
! ... allocations:
!debug
      print *,'before maxl'
      maxl=maxval(lpar)
      write(*,'(/a,i10)') 'max_Lpar = ', maxl
      print *,'maxl = ',maxl

      print *,'mom,nlsp,ne:',mom,nlsp,ne

      Allocate(fom(mom,nlsp,ne), fl(0:maxl), iop(nlsp,ne,3), om(mom), tmatr(mom), tmati(mom) )

      STOP

      fom = 0.d0;  iop = 0;  jop =0; om_top = 0.d0; om_sum = 0.d0;
      mem = mem + 2*mom*nlsp*ne + 2*maxl + nlsp*ne*3 + 6*mom

      mtr = np*(np+1)/2; if(ion.ne.0) mtr=np*(np-1)/2

      if(np.lt.ntarg) mtr = mtr + (ntarg-np)*ni

      Allocate( fail(mtr), coefa(mtr), coefb(mtr), met(mtr), ic(mtr), jc(mtr) )
      fail = 0.d0; coefa = 0.d0; coefb = 0.d0; met = -2;
      mem =  mem + 18*mtr

      Allocate(om_top(mtr,ne), om_sum(mtr,ne) )
      om_top = 0.d0; om_sum = 0.d0;
      mem =  mem + 4*mtr*ne

      Do itr1 = 1,np
      Do itr2 = itr1,ntarg
       itr = Index_TR(ion,itr1,itr2,np,ni)
       if(itr.eq.0) Cycle
       ic(itr) = itr1
       jc(itr) = itr2
      End do; End do

      write(*,'(/a,i5,a)') 'mtr   = ', mtr,' - number of transitions'

      maxll = maxval(lch)
      write(*,'(/a,i3,a)') 'max l = ', maxll

      Allocate(fll(0:maxll,2))
      mem =  mem + 4*maxll

      S = mem / (256.0*1024.0)
      write(*,'(/a,f10.2,a)') 'Memory required:  ', S, ' Mb'

      if(S.gt.50000.d0) Stop ' > 50 GB'

!-----------------------------------------------------------------------------------
! ... check continuation:

      if(Icheck_file(ccc).gt.0) then
       Open(nuc,file=ccc)
       read(nuc,*) mm;  if(mtr.ne.mm) Stop 'mtr <> mm'
       Do j=1,mtr
        read(nuc,*) ic(j),jc(j), fail(j), coefa(j), coefa(j), met(j)
       End do
       Close(nuc)
      end if

!----------------------------------------------------------------------
! ... read t-matrix and convert to channel omega:

      rewind(nut)

    3 if(form.eq.'a') then
       read(nut,*,end=4) e1,nopen,nom,ilsp; i1=np; i2=ni; kopen=nopen
       read(nut,*) (tmatr(i),tmati(i),i=1,nom)
      else
       read(nut,*,end=4) e1,nopen,kopen,ilsp,i1,i2,nj
       ntr = kopen*(kopen+1)/2
       read(nut,*) (tmatr(i),tmati(i),i=1,ntr)
       ktr = (nopen-kopen)*nj
       if(ktr.gt.0) read(nut,*) (tmatr(i),tmati(i),i=ntr+1,ntr+ktr)
      end if

      ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do
      if(ie.eq.0) go to 3

      iop(ilsp,ie,1)=nopen
      iop(ilsp,ie,2)=kopen
      iop(ilsp,ie,3)=nj

! ... partial omega:

      g = (2*lpar(ilsp)+1) * iabs(ispar(ilsp)) / 2.d0
      if(ispar(ilsp).eq.0) g = (lpar(ilsp)+1)/2.d0

      om=0.d0
      Do i=1,kopen
       Do j=1,i
        ij=(i-1)*i/2+j
        om(ij) = (tmatr(ij)*tmatr(ij)+tmati(ij)*tmati(ij))*g
       End do
      End do

      nom = kopen*(kopen+1)/2

      if(ktr.gt.0) then

       Do i=kopen+1,nopen
        Do j=1,nj;
         ij=nom+(i-kopen-1)*nj+j
         om(ij) = (tmatr(ij)*tmatr(ij)+tmati(ij)*tmati(ij))*g
        End do
       End do
       nom = nom + nj * (nopen-kopen)

      end if

      Do i=1,nom; fom(i,ilsp,ie) = om(i); End do

      go to 3
    4 Close(nut)

      Call CPU_time(t2)
      write(*,'(/a,f10.1,a)') 'read t-matrix = ',(t2-t1)/60,' min'


!-----------------------------------------------------------------------
! ... ouput files:

      if(len_trim(label).gt.0) top = trim(top)//'_'//trim(label)
      Call Read_aarg('top',top)
      Open(nuo,file=top)

      if(len_trim(label).gt.0) oms = trim(oms)//'_'//trim(label)
      Call Read_aarg('oms',oms)
      Open(nuq,file=oms)

      if(len_trim(label).gt.0) bad = trim(bad)//'_'//trim(label)
      Open(nub,file=bad)

!-----------------------------------------------------------------------
! ... cycle over energies:

      nbad = 0
      Do  ie=1,ne
       ek = e(ie)

!-----------------------------------------------------------------------
! ... cycle over transitions:

      Do itr1 = 1,ntarg
      Do itr2 = itr1,ntarg

       if(jtr1.ne.0.and.jtr1.ne.itr1) Cycle
       if(jtr2.ne.0.and.jtr2.ne.itr2) Cycle

       if(ek1.gt.0.d0.and.e1.lt.ek1) go to 1
       if(ek2.gt.0.d0.and.e1.gt.ek2) go to 1

       if(ek.lt.etarg(itr2)) Cycle

       itr = Index_TR(ion,itr1,itr2,np,ni)
       if(itr.eq.0) Cycle

       ! ... find met: type of transition

       if(IStarg(itr1).ne.IStarg(itr2)) then
        met(itr) = -1;  AS = 'exchange'                        ! exchange
       elseif(ISTARG(itr1).ne.0.and.   &
              iptarg(itr1).ne.iptarg(itr2).and. &
              ltarg(itr1)+ltarg(itr2).ne.0.and. &
              iabs(ltarg(itr1)-ltarg(itr2)).le.1) then
        met(itr) =  0;  AS = 'dipole'                          ! dipole, LS
       elseif(ISTARG(itr1).eq.0.and.   &
              iptarg(itr1).ne.iptarg(itr2).and. &
              ltarg(itr1)+ltarg(itr2).ne.0.and. &
              iabs(ltarg(itr1)-ltarg(itr2)).le.2) then
        met(itr) =  0;  AS = 'dipole'                          ! dipole, JK
       else
        met(itr) =  1;  AS = 'non-dipole'                      ! non-dipole
       end if

       write(pri,'(85(''-''))')
       write(pri,'(a,2i5,f15.8,5x,a)') 'transition ',itr1,itr2,ek,trim(AS)
       write(pri,'(85(''-''))')

!------------------------------------------------------------------------------------------

       e1 = ek-etarg(itr1);  ekap1=e1/zion**2;  kap1=sqrt(ekap1)
       e2 = ek-etarg(itr2);  ekap2=e2/zion**2;  kap2=sqrt(ekap2)

       if(IStarg(itr1).ne.0) then
        ll1 = lpar(nlsp)-ltarg(itr1)
        ll2 = lpar(nlsp)-ltarg(itr2)
        ll_max = min(ll1,ll2)
       else
        ll1 = NINT( real(jpar(nlsp)) - real(jtarg(itr1)) - 1) / 2
        ll2 = NINT( real(jpar(nlsp)) - real(jtarg(itr2)) - 1) / 2
        ll_max = min(ll1,ll2)
       end if
       lamda = ll_max - 1

       if(met(itr).eq.-1)  Call TOP_exchange
       if(met(itr).eq. 1)  Call TOP_exchange   ! TOP_nondipole        ???

       if(met(itr).eq.0 .and. method.eq.1 ) Call TOP_dipole1
       if(met(itr).eq.0 .and. method.eq.2 ) Call TOP_dipole2

!------------------------------------------------------------------------------------------

      End do         !  over  itr2
      End do         !  over  itr1

      End do         !  over ie

!-----------------------------------------------------------------------
!... fail information:

      open(nuc,file=ccc)
      rewind(nuc)
      write(nuc,'(i10,a)') mtr, ' - number of transitions'

      i = 0
      Do j= 1,mtr
       if(fail(j).eq.0.d0) then
        write(nuc,'(2i8,f16.8,2E15.5,2i5)') ic(j),jc(j), fail(j), coefa(j), coefb(j), met(j)
       else
        i = i + 1
        write(nuc,'(2i8,f16.8,2E15.5,2i5)') ic(j),jc(j), fail(j), coefa(j),coefb(j),met(j),i
       end if
      End do

      write(nuc,*) 'failed transitions: ',i
      Close(nuc)

      write(*,*) 'failed transitions: ',i

!----------------------------------------------------------------------
! ... output new 'topped' om:

      Do ie=1,ne

      mtr = np*(np+1)/2; if(ion.ne.0) mtr=np*(np-1)/2

      if(np.lt.ntarg) mtr = mtr + (ntarg-np)*ni

       nopen = Iopen(ntarg,e(ie),etarg)
       if(ion.gt.0.and.nopen.le.1) Cycle

       if(nopen.le.np) then
        if(ion.eq.0) nom = nopen*(nopen+1)/2
        if(ion.gt.0) nom = nopen*(nopen-1)/2
       else
        if(ion.eq.0) nom = np*(np+1)/2
        if(ion.gt.0) nom = np*(np-1)/2
        nom = nom + (nopen-np)*ni
       end if
       kopen = min(np,nopen)

       write(nuo,'(F10.6,5i8,a)')  e(ie),nom,nopen,kopen,np,ni, &
                              '   e(ie),nom,iopen,jopen,np,ni'
       write(nuo,'(5D16.8)') (om_top(i,ie),i=1,nom)

       write(nuq,'(F10.6,5i8,a)')  e(ie),nom,nopen,kopen,np,ni, &
                              '   e(ie),nom,iopen,jopen,np,ni'
       write(nuq,'(5D16.8)') (om_sum(i,ie),i=1,nom)

      End do        ! over ie

!----------------------------------------------------------------------

      write(*,*) 'failed energies: ',nbad

      Call CPU_time(t2)

      write(  *,'(/a,f10.1,a)') 'time = ',(t2-t1)/60,' min'
      write(pri,'(/a,f10.1,a)') 'time = ',(t2-t1)/60,' min'

 CONTAINS

!======================================================================
      Subroutine Extra1
!======================================================================
! ... extrapolation of dipole transitions
!----------------------------------------------------------------------

      if(coefa(itr).eq.0.d0) then
        je = ie-1
        x1 = e(je-2)
        y1 = om_top(itr,je-2)/(e(je-2)-etarg(i1))
        x2 = E(je)
        y2 = om_top(itr,je)/(e(je)-etarg(i1))
        y1 = y1/ log(x1)
        y2 = y2/ log(x2)
        b = (y1*x1-y2*x2)/(y2-y1)
        a = y1*(x1+b)
        coefa(itr) = a
        coefb(itr) = b
       end if

      x = e(ie);  a=coefa(itr); b=coefb(itr)
      om_top(itr,ie) = a*log(x)/(x+b)*(e(ie)-etarg(i1))

      End  Subroutine Extra1


!======================================================================
      Subroutine Extra2
!======================================================================
! ... extrapolation of non-dipole transitions
!----------------------------------------------------------------------

      if(coefa(itr).eq.0.d0) then
        i1 = max(1,ie-3)
        i2 = ie-1; if(ie.le.1) i2=ie
        coefa(itr) = sum(om_top(itr,i1:i2))/(i2-i1+1)
        coefb(itr) = 0.d0
      end if
      om_top(itr,ie) = coefa(itr)

      End  Subroutine Extra2

!======================================================================
      Subroutine Extra2_old
!======================================================================
! ... extrapolation of non-dipole transitions
!----------------------------------------------------------------------

      if(coefa(itr).eq.0.d0) then
        je = ie-1
        x1 = E(je-2)
        y1 = om_top(itr,je-2)/(e(je-2)-etarg(i1))
        x2 = E(je)
        y2 = om_top(itr,je)/(e(je)-etarg(i1))
        b = (y1*x1-y2*x2)/(y1-y2)
        a = y1*(x1-b)
        coefa(itr) = a
        coefb(itr) = b
      end if

      x = e(ie);  a=coefa(itr);  b=coefb(itr)
      om_top(itr,ie) = a/(x-b)*(e(ie)-etarg(i1))

      End  Subroutine Extra2_old



!======================================================================
      Subroutine TOP_exchange
!======================================================================
! ... top-up for exchange transitions, i.e., no  top-up
!----------------------------------------------------------------------
      Call Get_fl

      S=SUM(fl);  om_sum(itr,ie)=S;  om_top(itr,ie)=S

      f1=0.d0; f2=0.d0; f3=0.d0
      Do il=0,maxl
       if(fl(il).eq.0.d0) Cycle
       f1=f2; f2=f3; f3=fl(il); C = f2/f3; if(C.lt.1) C=1.d0
       if(jtr1.ne.0) write(pri,'(i2,a1,E15.5,f10.3)') il,'.',f3,C-1.d0
      End do

      if(jtr1.ne.0) write(pri,'(18(''-''))')
      write(pri,'(a3,E15.5)') 'SUM',S
      write(pri,'(18(''-''))')
      write(pri,'(a3,E15.5,a)') 'TOP',S, ' -  exchange, no top-up'

      End Subroutine TOP_exchange



!======================================================================
      Subroutine TOP_nondipole
!======================================================================
! ... top-up for exchange transitions, geometric series
!----------------------------------------------------------------------
      Call Get_fl

      S=SUM(fl);  om_sum(itr,ie)=S;  om_top(itr,ie)=S

      f1=0.d0; f2=0.d0; f3=0.d0
      Do il=0,maxl
       if(fl(il).eq.0.d0) Cycle
       f1=f2; f2=f3; f3=fl(il); C = f2/f3; if(C.lt.1) C=1.d0; if(C.gt.10) C=1.d0
       if(jtr1.ne.0) write(pri,'(i2,a1,E15.5,f10.3)') il,'.',f3,C-1.d0
      End do

      if(jtr1.ne.0) write(pri,'(18(''-''))')
      write(pri,'(a3,E15.5)') 'SUM',S
      write(pri,'(18(''-''))')

      if(fail(itr).gt.0.d0.and.ek.ge.fail(itr)) then
       Call Extra2
       write(pri,'(a3,E15.5,a)') 'TOP',om_top(itr,ie), ' -  extrapolated'
      elseif(S*eps_tail.gt.f1+f2+f3) then
       write(pri,'(a3,E15.5,a)') 'TOP',S, ' -  small correction, no top-up'
      else

       x1=f1/f2-1.d0; x2=f2/f3-1.d0; xx=(x1+x2)/2
       xe=(EK-Etarg(itr1))/(EK-Etarg(itr2)) - 1.d0

       if(fail(itr).eq.0.and.xx.lt.eps_x)  then
        if(ek.gt.ekk) then
         fail(itr)=ek
        else
         write(nub,'(f10.6,2i5)') ek,itr1,itr2; nbad=nbad+1; Return
        end if
       end if

       if(fail(itr).eq.0) then
        x = xx;   S = S + f3/x;   om_top(itr,ie)=S
        write(pri,'(a3,E15.5,3F10.3,a)') 'TOP',S,xx,xe,x,'  average, energy-drivenv, chosen '
       else
        Call Extra2
        write(pri,'(a3,E15.5,a)') 'TOP',om_top(itr,ie),' -  extrapolated'
       end if

      end if

      End Subroutine TOP_nondipole


!======================================================================
      Subroutine TOP_dipole1
!======================================================================
! ... top-up for exchange transitions, geometric series
!----------------------------------------------------------------------
      Call Get_fl

      kbad=0
      Do ilsp = 1,nlsp;  if(fl(lpar(ilsp)).eq.0.d0) kbad=1;  End do

      if(kbad.gt.0) then
       write(nub,'(f10.6,2i5)') ek,itr1,itr2; nbad=nbad+1
      end if

      S=SUM(fl);  om_sum(itr,ie)=S;  om_top(itr,ie)=S

      f1=0.d0; f2=0.d0; f3=0.d0
      Do il=0,maxl
       if(fl(il).eq.0.d0) Cycle
       f1=f2; f2=f3; f3=fl(il); C = f2/f3; if(C.lt.1.d0) C=1.d0
       if(jtr1.ne.0) write(pri,'(i2,a1,E15.5,f12.3)') il,'.',f3,C-1.d0
      End do

      if(jtr1.ne.0) write(pri,'(18(''-''))')
      write(pri,'(a3,E15.5)') 'SUM',S
      write(pri,'(18(''-''))')

      if(fail(itr).gt.0.d0.and.ek.ge.fail(itr)) then

       Call Extra1
       write(pri,'(a3,E15.5,a)') 'TOP',om_top(itr,ie), ' -  extrapolated'

      elseif(S*eps_tail.gt.f1+f2+f3) then

       write(pri,'(a3,E15.5,a)') 'TOP',S, ' -  small correction, no top-up'

      else

       x1=f1/f2-1.d0; x2=f2/f3-1.d0; xx=(x1+x2)/2
       xe=(EK-Etarg(itr1))/(EK-Etarg(itr2)) - 1.d0

       if(fail(itr).eq.0.and.xx.lt.eps_x.and.ie.gt.3)  then
        if(ek.gt.ekk) then
         fail(itr)=ek
        else
         write(nub,'(f10.6,2i5)') ek,itr1,itr2; nbad=nbad+1
        end if
       end if

       if(fail(itr).eq.0) then
        x = xx;   S = S + f3/x;   om_top(itr,ie)=S
        write(pri,'(a3,E15.5,3F10.3,a)') 'TOP',S,xx,xe,x,'  average, energy-drivenv, chosen '
       else
        Call Extra1
        write(pri,'(a3,E15.5,a)') 'TOP',om_top(itr,ie),' -  extrapolated'
       end if

      end if

      End Subroutine TOP_dipole1



!======================================================================
      Subroutine TOP_dipole2
!======================================================================
! ... top-up for exchange transitions, CBE
!----------------------------------------------------------------------
      Call Get_fll

      S=SUM(fll);  om_sum(itr,ie)=S;  om_top(itr,ie)=S

      if(jtr1.ne.0) then
       Do l=1,ll_max+1
        eps = 1.d-3
!        F1 = Fdip0(ekap1,l,  ekap2,l-1,eps,ifail1)
        F1 = IZ0_lamda(1,kap1,l,kap2,l-1)
!        F2 = Fdip0(ekap1,l-1,ekap2,l,  eps,ifail2)
        F2 = IZ0_lamda(1,kap1,l-1,kap2,l)
        om1 = 16.d0/3 * l * OS(itr1,itr2) * F1**2
        om2 = 16.d0/3 * l * OS(itr1,itr2) * F2**2
         write(pri,'(i2,a1,4E15.5)') l,'.',fll(l,:),om1,om2
       End do
       write(pri,'(18(''-''))')
      end if

      write(pri,'(a3,E15.5)') 'SUM',S
      write(pri,'(85(''-''))')

      S = SUM( fll(1:lamda,:) )  +  fll(lamda+1,2)

!      F1 = Fdip0(ekap1,lamda,  ekap2,lamda+1,eps,ifail1)
      F1 = IZ0_lamda(1,kap1,lamda,kap2,lamda+1)
!      F2 = Fdip0(ekap1,lamda+1,ekap2,lamda,  eps,ifail2)
      F2 = IZ0_lamda(1,kap1,lamda+1,kap2,lamda)
      om1 = 16.d0/3 * (lamda+1) * OS(itr1,itr2) * F1**2
      om2 = 16.d0/3 * (lamda+1) * OS(itr1,itr2) * F2**2
      C = (om2-om1) * (real(ion**2)/(lamda+1) + e1) / (e1-e2)

!      if(ifail1.eq.0.and.ifail2.eq.0) then
       S = S + C;  om_top(itr,ie) = S
       write(pri,'(a3,E15.5,a,i3,E15.5)') &
          'TOP',S,' - CBE approximation, lamda = ',lamda, C
!      else
!        if(ek.gt.ekk) then
!         fail(itr)=ek
!         S = S + C; om_top(itr,ie) = S
!         write(pri,'(a3,E15.5,a,i3)') 'TOP',S,' - CBE approximation failed'
!        else
!         write(nub,'(f10.6,2i5)') ek,itr1,itr2; nbad=nbad+1; Return
!        end if
!      end if


      End Subroutine TOP_dipole2


!======================================================================
      Subroutine TOP_dipole3
!======================================================================
! ... top-up for exchange transitions, CBE ???
!----------------------------------------------------------------------
      Call Get_fll

      S=SUM(fll);  om_sum(itr,ie)=S;  om_top(itr,ie)=S

      Do l=1,ll_max

       eps = 1.d-5
       F1 = Fdip0(e1,l,e2,l-1,eps,ifail)
       F2 = Fdip0(e1,l-1,e2,l,eps,ifail)
       if(ifail.ne.0) write(*,'(a,i5,2f10.5)') 'Fdip0 fail: ',l,e1,e2
       om1 = 16.d0/3 * l * OS(itr1,itr2) * F2**2
       om2 = 16.d0/3 * l * OS(itr1,itr2) * F1**2
       if(jtr1.ne.0) write(pri,'(i2,a1,4E15.5)') l,'.',fll(l,:),om1,om2

      End do
      if(jtr1.ne.0) write(pri,'(18(''-''))')
      write(pri,'(a3,E15.5)') 'SUM',S
      write(pri,'(85(''-''))')

      S = SUM(fll(1:ll_max,:))

      lamda = ll_max
      C = (1+lamda**2*e2)*om1 - (1+lamda**2*e1)*om2
      C = C /lamda**2 / (e1-e2)

      S = S + C;  om_top(itr,ie) = S

      write(pri,'(a3,E15.5,a,i3)') 'TOP',S,' - CBE approximation, lamda = ',lamda


      End Subroutine TOP_dipole3


!======================================================================
      Subroutine GET_fl
!======================================================================

      fl = 0.d0
      Do ilsp = 1,nlsp;  l = lpar(ilsp);  nopen = iop(ilsp,ie,1)
                                          kopen = iop(ilsp,ie,2)

        Do i=1,kopen;   if(iptar(ilsp,i).ne.itr2) Cycle
         Do j=1,i;      if(iptar(ilsp,j).ne.itr1) Cycle
          ij=(i-1)*i/2+j
          s = fom(ij,ilsp,ie)
          if(itr1.eq.itr2.and.i.ne.j) s = s + s
          fl(l) = fl(l) +  s
         End do
        End do

        if(nopen.gt.kopen) then

         nj = iop(ilsp,ie,3)
         ntr = kopen*(kopen+1)/2
         Do i=kopen+1,nopen;  if(iptar(ilsp,i).ne.itr2) Cycle
          Do j=1,nj;          if(iptar(ilsp,j).ne.itr1) Cycle
           ij=ntr+(i-kopen-1)*nj+j
           s = fom(ij,ilsp,ie)
           if(itr1.eq.itr2.and.i.ne.j) s = s + s
           fl(l) = fl(l) +  s
          End do
         End do

        end if

      End do  ! over ilsp

      End Subroutine GET_fl


!======================================================================
      Subroutine GET_fll
!======================================================================

      fll = 0.d0
      Do ilsp = 1,nlsp;  nopen = iop(ilsp,ie,1);  kopen = iop(ilsp,ie,2)

        Do i=1,kopen;   if(iptar(ilsp,i).ne.itr2) Cycle; l1=lch(ilsp,i)
         Do j=1,i;      if(iptar(ilsp,j).ne.itr1) Cycle; l2=lch(ilsp,j)
          ij=(i-1)*i/2+j
          s = fom(ij,ilsp,ie)
          if(itr1.eq.itr2.and.i.ne.j) s = s + s
          if(l1.gt.l2) fll(l1,1) = fll(l1,1) + s
          if(l2.gt.l1) fll(l2,2) = fll(l2,2) + s
         End do
        End do

        if(nopen.gt.kopen) then

         nj = iop(ilsp,ie,3)
         ntr = kopen*(kopen+1)/2
         Do i=kopen+1,nopen;  if(iptar(ilsp,i).ne.itr2) Cycle; l1=lch(ilsp,i)
          Do j=1,nj;          if(iptar(ilsp,j).ne.itr1) Cycle; l2=lch(ilsp,j)
           ij=ntr+(i-kopen-1)*nj+j
           s = fom(ij,ilsp,ie)
           if(itr1.eq.itr2.and.i.ne.j) s = s + s
           if(l1.gt.l2) fll(l1,1) = fll(l1,1) + s
           if(l2.gt.l1) fll(l2,2) = fll(l2,2) + s
          End do
         End do

        end if

      End do  ! over ilsp

      End Subroutine GET_fll

      End   !  program sec_top


!======================================================================
      Subroutine inf_top_omb
!======================================================================
!     provide screen information about sec_top_TM utility
!----------------------------------------------------------------------
      Character(80) :: A

      Call get_command_argument(1,A)
      if(A.ne.'?') Return

      write(*,'(a)') &
'                                                                  ',&
'Top up of the collision strengths based on the geometric series   ',&
'approximation (method 1) and CBE approximation (method 2) for     ',&
'dipole transitions.                                               ',&
'                                                                  ',&
'     zarm.tma (zarm.tmb) --> zarm.tmb_top                         ',&
'                                                                  ',&
'Optional arguments:                                               ',&
'                                                                  ',&
'           par  - file with partial omega"s    [zarm.omb_par]     ',&
'           top  - file with topped up omega"s  [zarm.omb_top]     ',&
'           ek1, ek2 - energy interval in Ry                       ',&
'           tail - tolerence for top-up contribution  [0.001]      ',&
'           x    - tolerence for geom. series parameter [0.05]     ',&
'           method [1] - if =2, CBE is used for dipole transitions ',&
'                                                                  '
      Stop ' '

      End Subroutine inf_top_omb



