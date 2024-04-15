!======================================================================
!     UTILITY      S E C _ T O P
!
!     zarm.omb_par --> zarm.omb_top
!
!======================================================================
!     generate top-up omegas for all included transitions
!     (based on the information in 'zarm.omb_par')
!
!     Call as:  sec_top  [par= top= ek1= ek2= eps_tail= eps_fail= eps_x=]
!---------------------------------------------------------------------

      Use target; Use channels

      Implicit real(8) (a-h,o-z)

      Real(8), Allocatable ::  fl(:), om(:), e(:),fail(:),coefa(:),coefb(:)
      Real(8), Allocatable ::  om_sum(:,:), om_top(:,:)
      Real(8), Allocatable ::  fom(:,:,:)

      Integer, Allocatable ::  iop(:),jop(:), met(:), ic(:), jc(:)

      Integer :: ke = 20000  !  initial number of energies

      Real(8) :: eps_tail = 0.001
      Real(8) :: eps_x    = 0.025

      Integer :: np=0, ni=0
!----------------------------------------------------------------------
! ... files:
      Character(80) :: AF, label=' '
      Integer :: nup=11;  Character(80) :: targ  = 'target'
      Integer :: nut=12;  Character(80) :: par   = 'zarm.omb_par'
      Integer :: nuo=14;  Character(80) :: top   = 'zarm.omb_top'
      Integer :: nuq=15;  Character(80) :: oms   = 'zarm.omb'
      Integer :: pri=16;  Character(80) :: out   = 'sec_top_omb.log'
      Integer :: nuc=17;  Character(80) :: ccc   = 'sec_top_coef_fail'
      Integer :: nub=18;  Character(80) :: bad   = 'bad_energies'

      Integer :: nua=99;  ! scratch file

      Call CPU_time(t1)

      Call Inf_top_omb

      Open(pri,file=out)

!----------------------------------------------------------------------
! ... target information:

      Call Check_file(targ)
      Open(nup,file=targ)
      Call R_target(nup)
      Call R_channels(nup)
      Close(nup)
      e1=etarg(1); Do i=1,ntarg; etarg(i)=(etarg(i)-e1)*2; End do
      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion

      Z = nz
      AWT = 0.d0;  Call Read_rarg('AWT',AWT)
      Call Conv_au (Z,AWT,au_cm,au_eV,pri)
      Ry = au_eV/2.d0

!-----------------------------------------------------------------------
! ... define other parameters:

      ek1 = 0.d0; Call Read_rarg('ek1',ek1)
      ek2 = 0.d0; Call Read_rarg('ek2',ek2)
      ekk = 0.d0; Call Read_rarg('ekk',ekk)

      eps_tail=0.001; Call Read_rarg('tail',eps_tail)
      eps_x=0.025;    Call Read_rarg('x',eps_x)

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

!----------------------------------------------------------------------
! ... find energies:

      if(len_trim(label).gt.0) par = trim(par)//'_'//trim(label)
      Call Read_aarg('par',par)
      Call Check_file(par)
      Open(nut,file=par)

      me = ke; Allocate(e(me))

      ne=0; mom=0
    1 read(nut,*,end=2) e1,ilsp,nom,k1,k2,i1,i2
      read(nut,*) (S,i=1,nom)

      if(ek1.gt.0.d0.and.e1.lt.ek1) go to 1
      if(ek2.gt.0.d0.and.e1.gt.ek2) go to 1

      if(np.eq.0) np = i1
      if(ni.eq.0) ni = i2
      if(np.ne.i1) Stop 'diferent np'
      if(ni.ne.i2) Stop 'diferent ni'

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
    2 write(*,'(/a,i5,a)') 'ne =',ne,' - number of energies'
      if(ne.eq.0) Stop 'nothing to do !'

      Call Rsort(ne,e)

      write(*,'(/a,2f15.6,a)') 'e(1),e(ne)', e(1), e(ne), ' - energy interval'
      write(*,'(/a,i5,a)') 'mom = ',mom,' - maximum matrix dimension'

      write(*,'(/a,i5,a)') 'np =',np,' - number of physical states'
      write(*,'(/a,i5,a)') 'ni =',ni,' - number of "ionization" states'

!----------------------------------------------------------------------
! ... allocations:

      maxl=maxval(lpar)
      write(*,'(/a,3i10/)') 'maxl = ', maxl

      S = 8.0 * mom * (nlsp+3)  / (1024.0*1024.0);  S = S * ne
      write(*,'(/a,f10.2,a/)') 'Memory required:  ', S, ' Mb'
      if(S.gt.50000.d0) Stop ' > 50 GB'

      Allocate(fom(mom,nlsp,ne), fl(0:maxl), iop(ne), jop(ne), om(mom) )

      Allocate(om_top(mom,ne), om_sum(mom,ne) )

      fom = 0.d0;  iop = 0;  jop =0; om_top = 0.d0; om_sum = 0.d0;

      mmm = np*(np+1)/2; if(ion.ne.0) mmm=np*(np-1)/2

      if(np.lt.ntarg) mmm = mmm + (ntarg-np)*ni

      Allocate( fail(mmm), coefa(mmm), coefb(mmm), met(mmm), ic(mmm), jc(mmm) )
      fail = 0.d0; coefa = 0.d0; coefb = 0.d0; met = -2;

      Do itr1 = 1,np
      Do itr2 = itr1,ntarg
       itr = Index_TR(ion,itr1,itr2,np,ni)
       if(itr.eq.0) Cycle
       ic(itr) = itr1
       jc(itr) = itr2
      End do; End do

!-----------------------------------------------------------------------------------
! ... check continuation:

      if(Icheck_file(ccc).gt.0) then
       Open(nuc,file=ccc)
       read(nuc,*) mm;  if(mmm.ne.mm) Stop 'mmm <> mm'
       Do j=1,mmm
        read(nuc,*) ic(j),jc(j), fail(j), coefa(j), coefa(j), met(j)
       End do
       Close(nuc)
      end if

!----------------------------------------------------------------------
! ... read partial OM:

      rewind(nut)
    3 read(nut,*,end=4) e1,ilsp,nom,k1,k2,i1,i2
      read(nut,*) (om(i),i=1,nom)

      ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do
      if(ie.eq.0) go to 3

      if(iop(ie).eq.0) iop(ie)=k1
      if(jop(ie).eq.0) jop(ie)=k2

      if(iop(ie).ne.k1) then
       write(*,'(i5,f14.6,10i5)') ilsp, e1, ie, k1, iop(ie), IOPEN(ntarg,e1,etarg)
       Stop 'different iop'
!       iop(ie)=-10
      end if

      if(jop(ie).ne.k2) Stop 'different jop'

      Do i=1,nom; fom(i,ilsp,ie) = om(i); End do

      go to 3
    4 Close(nut)

      Call CPU_time(t2)
      write(*,'(a,f10.1,a)') 'read omega = ',(t2-t1)/60,' min'
      write(pri,'(a,f10.1,a)') 'read omega = ',(t2-t1)/60,' min'

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
       om = 0.d0

!       if(iop(ie).lt.0) Cycle

!-----------------------------------------------------------------------
! ... cycle over transitions:

      ntr1 = iop(ie);  ntr = ntr1*(ntr1+1)/2
      if(ion.ne.0)     ntr = ntr1*(ntr1-1)/2
      ntr2 = jop(ie);  nom = ntr + (ntr2-ntr1)*ni

      Do itr1 = 1,ntr1
      Do itr2 = itr1,ntr2

       if(jtr1.ne.0.and.jtr1.ne.itr1) Cycle
       if(jtr2.ne.0.and.jtr2.ne.itr2) Cycle

       if(ek.lt.etarg(itr2)) Cycle

       itr = Index_TR(ion,itr1,itr2,np,ni)
       if(itr.eq.0) Cycle

       ! ... find met:

       if(IStarg(itr1).ne.IStarg(itr2)) then
        met(itr) = -1                                          ! exchange
       elseif(ISTARG(itr1).ne.0.and.   &
              iptarg(itr1).ne.iptarg(itr2).and. &
              ltarg(itr1)+ltarg(itr2).ne.0.and. &
              iabs(ltarg(itr1)-ltarg(itr2)).le.1) then
        met(itr) =  0                                          ! dipole, LS
       elseif(ISTARG(itr1).eq.0.and.   &
              iptarg(itr1).ne.iptarg(itr2).and. &
              ltarg(itr1)+ltarg(itr2).ne.0.and. &
              iabs(ltarg(itr1)-ltarg(itr2)).le.2) then
        met(itr) =  0                                          ! dipole, JK
       else
        met(itr) =  1                                          ! non-dipole
       end if

       write(pri,'(15(''-''))')
       write(pri,'(a,3i5,f15.8,i5)') 'transition ',itr,itr1,itr2,ek,met(itr)
       write(pri,'(15(''-''))')

       fl = 0.d0
       Do ilsp = 1,nlsp;  l = lpar(ilsp);  fl(l) = fl(l) + fom(itr,ilsp,ie);  End do

       kbad=0
       Do ilsp = 1,nlsp;  if(fl(lpar(ilsp)).eq.0.d0.and.ek.gt.0.23) kbad=1;  End do
       if(kbad.gt.0) then; write(nub,'(f10.6,2i5)') ek,itr1,itr2; nbad=nbad+1; Cycle;  end if

       S=SUM(fl);  om_sum(itr,ie)=S;  om_top(itr,ie)=S;  if(S.eq.0.d0) Cycle

       f1=0.d0; f2=0.d0; f3=0.d0
       Do il=0,maxl
        if(fl(il).eq.0.d0) Cycle
        f1=f2; f2=f3; f3=fl(il)
        if(f2.eq.0.d0) f2=f3
        if(jtr1.ne.0) write(pri,'(i2,a1,D12.4,f10.3)') il,'.',f3,f2/f3-1.d0
       End do
       write(pri,'(15(''-''))')
       write(pri,'(a3,D12.4)') 'SUM',S
       write(pri,'(15(''-''))')


       if(IStarg(itr1).ne.IStarg(itr2)) then

        write(pri,'(a3,D12.4,a)') 'TOP',S, ' -  exchange, no top-up'

       elseif(fail(itr).gt.0.d0.and.ek.ge.fail(itr)) then

        if(met(itr).eq.0) Call Extra1
        if(met(itr).eq.1) Call Extra2

        write(pri,'(a3,D12.4,a)') 'TOP',om_top(itr,ie), ' -  extrapolated'

       elseif(S*eps_tail.gt.f1+f2+f3) then

        write(pri,'(a3,D12.4,a)') 'TOP',S, ' -  small correction, no top-up'

       else

        x1=f1/f2-1.d0; x2=f2/f3-1.d0; xx=(x1+x2)/2
        xe=(EK-Etarg(itr1))/(EK-Etarg(itr2)) - 1.d0

        if(fail(itr).eq.0.and.xx.lt.eps_x)  then
         if(ek.gt.ekk) then
          fail(itr)=ek
         else
          write(nub,'(f10.6,2i5)') ek,itr1,itr2; nbad=nbad+1; Cycle;
         end if
        end if

        if(fail(itr).eq.0) then
         x = xx
         S = S + f3/x;   om_top(itr,ie)=S
         write(pri,'(a3,D12.4,3F10.3,a)') 'TOP',S,xx,xe,x,'  average, energy-drivenv, chosen '
        else
         if(met(itr).eq.0) Call Extra1
         if(met(itr).eq.1) Call Extra2
         write(pri,'(a3,D12.4,a)') 'TOP',om_top(itr,ie),' -  extrapolated'
        end if

       end if

       write(pri,'(15(''-''))')

      End do         !  over  itr2
      End do         !  over  itr1

      End do        ! over ie

!-----------------------------------------------------------------------
!... fail information:

      open(nuc,file=ccc)
      rewind(nuc)
      write(nuc,'(i10,a)') mmm, ' - number of transitions'

      i = 0
      Do j= 1,mmm
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
       if(iop(ie).lt.0) Cycle

       ntr1 = iop(ie);  ntr = ntr1*(ntr1+1)/2
       if(ion.ne.0)     ntr = ntr1*(ntr1-1)/2
       ntr2 = jop(ie);  nom = ntr + (ntr2-ntr1)*ni

       write(nuo,'(F10.6,5i8,a)')  e(ie),nom,iop(ie),jop(ie),np,ni, &
                              '   e(ie),nom,iopen,jopen,np,ni'
       write(nuo,'(5D16.8)') (om_top(i,ie),i=1,nom)

       write(nuq,'(F10.6,5i8,a)')  e(ie),nom,iop(ie),jop(ie),np,ni, &
                              '   e(ie),nom,iopen,jopen,np,ni'
       write(nuq,'(5D16.8)') (om_sum(i,ie),i=1,nom)

      End do        ! over ie

!----------------------------------------------------------------------

      write(*,*) 'failed energies: ',nbad

      Call CPU_time(t2)

      write(*,'(a,f10.1,a)') 'time = ',(t2-t1)/60,' min'
      write(pri,'(a,f10.1,a)') 'time = ',(t2-t1)/60,' min'

 CONTAINS

!======================================================================
      Subroutine Extra1
!======================================================================
! ... extrapolation of dipole transitions
!----------------------------------------------------------------------

      Real(8) :: exp=2.718281828, z


      if(coefa(itr).eq.0.d0) then
        je = ie-1
        z  = e(je) /(etarg(itr2)-etarg(itr1))
        a  = om_top(itr,je)/ log(z)
write(*,'(2E12.5)') e(ie), om_top(itr,ie-1)
write(*,'(a,2E12.5)') 'z,a',z,a
        coefa(itr) = a
        coefb(itr) = 0.d0
      end if

      z  = e(ie) /(etarg(itr2)-etarg(itr1))
      a = coefa(itr)
      om_top(itr,ie) = a*log(z)

      End  Subroutine Extra1


!======================================================================
      Subroutine Extra2
!======================================================================
! ... extrapolation of non-dipole transitions
!----------------------------------------------------------------------
!     write(688,'(a)') 'wanted to call Extra2'
!     return
      if(coefa(itr).eq.0.d0) then
        je = ie-1
        if (je-2.le.0) then
         write (689,*) je-2
!       stop
        else
         x1 = E(je-2)
         y1 = om_top(itr,je-2)/(e(je-2)-etarg(itr1))
         x2 = E(je)
         y2 = om_top(itr,je)/(e(je)-etarg(itr1))
         b = (y1*x1-y2*x2)/(y1-y2)
         a = y1*(x1-b)
         coefa(itr) = a
         coefb(itr) = b
        endif
      end if

      x = e(ie);  a=coefa(itr);  b=coefb(itr)
      om_top(itr,ie) = a/(x-b)*(e(ie)-etarg(itr1))

      End  Subroutine Extra2

      End   !  program sec_top


!======================================================================
      Subroutine inf_top_omb
!======================================================================
!     provide screen information about sec_top_TM utility
!----------------------------------------------------------------------

      Character :: A=' '

      Call get_command_argument(1,A)
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                  ',&
'SEC_TOP tops up the zarm.omb_par file based on the                ',&
'geometric series approximation                                    ',&
'                                                                  ',&
'Arguments: par  - file with partial omega"s    [zarm.omb_par]     ',&
'           top  - file with topped up omega"s  [zarm.omb_top]     ',&
'           ek1, ek2 - energy interval in Ry                       ',&
'           eps_tail - tolerence for top-up contribution  [0.02]   ',&
'           eps_x    - tolerence for geom. series parameter [0.001]'
      Stop ' '

      End Subroutine inf_top_omb



