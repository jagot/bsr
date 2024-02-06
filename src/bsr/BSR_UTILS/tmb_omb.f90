!======================================================================
!     UTILITY  tmb_omb
!======================================================================
!
!     T-matrix  --> omega, omega_par
!
!     zarm.tmb  --> zarm.omb,  zarm.omb_par   (default names)
!
!     Call:  tmb_omb
!
!     This utility also provides list of energies (tmb_omb.log)
!     for which the collection of T-matrixes is incomplete or dublicated
!
!======================================================================

      Use target; Use channels

      Implicit real(8) (a-h,o-z)

      Parameter (ke=1000)

      Real(8), allocatable :: e(:), tmar(:),tmai(:), om(:)
      Real(8), allocatable :: fl(:,:)
      Integer, allocatable :: npl(:), ipe(:), iop(:), jop(:)
      Integer, allocatable :: ipl(:,:)

! ... files:

      Character(80) :: param = 'target';        Integer :: nup = 1
      Character(80) :: oname = 'zarm.omb';      Integer :: nuo = 2
      Character(80) :: qname = 'zarm.omb_par';  Integer :: nuq = 3
      Character(80) :: tname = 'zarm.tmb';      Integer :: nut = 4
      Character(80) :: dname = 'zarm.dmb';      Integer :: nud = 9    ! corrected data: without repeating

      Character(80) :: AF_log = 'tmb_omb.log';  Integer :: pri = 6
      Character(80) :: AF

      Integer :: nua=99   ! scratch file

!----------------------------------------------------------------------

      Call get_command_argument(1,AF)

      if(AF.eq.'?'.or. AF.eq.'!') then

        write(*,'(/a)') 'tmb_omb calculate OMEGA data '
        write(*,'(/a)') 'zarm.tmb + target  -->  zarm.omb,  zarm.omb_par'
        write(*,'(/a)') 'Call as:  tmb_omb [tm=.. om=.. par=..]'
        write(*,'(/a)') 'tmb_omb.log contains abcent or access informatiom'
        write(*,*)
        Stop ' '

      end if

!----------------------------------------------------------------------
! ... target information:

      Call Check_file(param)
      Open(nup,file=param);
      Call R_target(nup)
      Call R_channels(nup)
      np=ntarg; Call Read_ipar(nup,'np',np); Call Read_iarg('np',np)
      ni=ntarg; Call Read_ipar(nup,'ni',ni); Call Read_iarg('ni',ni)
      Close(nup)
      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion
      E1 = etarg(1); etarg = (etarg-E1)*2.d0

      Call Read_aarg('tm',tname)
      Call Check_file(tname)
      Call Read_aarg('om',oname)
      Call Read_aarg('par',qname)

      id=0; Call Read_iarg('d',id)
      if(id.eq.0) nud=0
      Open(nut,file=tname)
      if(nud.gt.0) Open(nud,file=dname)

      mdim=mch*(mch+1)/2; Allocate(tmar(mdim),tmai(mdim))

      mom=ntarg*(ntarg+1)/2; Allocate(om(mom))

!----------------------------------------------------------------------
! ... find energies:

      me = ke; Allocate(e(me))

      ne=0
    1 read(nut,*,end=2) e1,nopen,kopen,ilsp,i1,i2,nj
      read(nut,*) ((S1,S2,j=1,i),i=1,kopen)
      if(nopen.gt.kopen.and.nj.gt.0) read(nut,*) &
                  ((S1,S2,j=1,nj),i=kopen+1,nopen)
      if(np.eq.0) np = i1
      if(ni.eq.0) ni = i2
      if(np.ne.i1) Stop 'diferent np'
      if(ni.ne.i2) Stop 'diferent ni'

      ie=0;  Do i=1,ne; if(e1.ne.e(i)) Cycle; ie=i; Exit; End do

      if(ie.eq.0) then; ne=ne+1; e(ne)=e1;  end if

      if(ne.eq.me) then
       open(nua,form='UNFORMATTED',status='SCRATCH')
       rewind(nua); write(nua) (e(i),i=1,ne)
       Deallocate(e); me=me+ke; Allocate(e(me))
       rewind(nua); read(nua) (e(i),i=1,ne)
      end if

      go to 1
    2 write(*,*) ' ne =',ne
      if(ne.eq.0) Stop ' '


! ... define dimensions (iop,jop):

      Allocate(iop(ne),jop(ne),npl(ne)); iop=0; jop=0; npl=0
      mom=0
      Do ie=1,ne
       i = Iopen(ntarg,e(ie),ETARG)
       iop(ie) = i; jop(ie)=i
       if(i.gt.np) iop(ie)=np
       npl(ie) = iop(ie)*(iop(ie)+1)/2
       if(jop(ie).gt.iop(ie)) npl(ie)=npl(ie) + ni*(jop(ie)-iop(ie))
       if(npl(ie).gt.mom) mom=npl(ie)
      End do

      Allocate(fl(mom,ne),ipl(nlsp,ne)); fl=0.d0; ipl=0

!----------------------------------------------------------------------
! ... read and transform T-matrices:

      Open(nuo,file=oname)
      Open(nuq,file=qname)

      rewind(nut); if(nud.gt.0) rewind(nud)
    3 read(nut,*,end=4) ee,nopen,kopen,ilsp,i1,i2,nj

      ntr = kopen*(kopen+1)/2
      read(nut,*) (tmar(i),tmai(i),i=1,ntr)
      ktr = (nopen-kopen)*nj
      if(ktr.gt.0) read(nut,*) (tmar(i),tmai(i),i=ntr+1,ntr+ktr)

      if(kopen.gt.nopen)     Stop ' kopen > nopen'
      if(ilsp.gt.nlsp)       Stop ' ilsp  > nlsp'
      if(ilsp.lt.0)          Stop ' ilsp < 0'
      if(nopen.gt.nch(ilsp)) then
         write(*,*) ' nopen > nch(ilsp)', nopen,nch(ilsp), ilsp,ee
         go to 3
      end if

      ! Stop ' nopen > nch(ilsp)'

      ie=0
      Do i=1,ne; if(ee.ne.e(i)) Cycle; ie=i; Exit; End do
      if(ie.eq.0) Stop 'unknown energy'

      ipl(ilsp,ie) = ipl(ilsp,ie) + 1
      if(ipl(ilsp,ie).gt.1) go to 3       !  ???

      if(nud.gt.0) then
       write(nud,'(F10.6,6i6,a)') ee,nopen,kopen,ilsp,i1,i2,nj, &
                              '   ee,nopen,kp,ilsp,np,ni,nj'
       write(nud,'(6d16.8)') (tmar(i),tmai(i),i=1,ntr)
       if(ktr.gt.0) write(nud,'(6d16.8)') (tmar(i),tmai(i),i=ntr+1,ntr+ktr)
      end if

!      if(iptar(ilsp,kopen).gt.iop(ie)) Stop 'kopen > iop(ne)'
!      if(iptar(ilsp,nopen).gt.jop(ie)) Stop 'nopen > jop(ne)'

      if(iptar(ilsp,kopen).gt.iop(ie)) write(*,*) 'kopen > iop(ne)', ee,nopen,kopen,ilsp,iop(ie)
      if(iptar(ilsp,nopen).gt.jop(ie)) write(*,*) 'nopen > jop(ne)', ee,nopen,kopen,ilsp,jop(ie)

! ... partial omega:

      g = (2*lpar(ilsp)+1) * iabs(ispar(ilsp)) / 2.d0
      if(ispar(ilsp).eq.0) g = (lpar(ilsp)+1)/2.d0

      om=0.d0
      Do i=1,kopen;    itr1=iptar(ilsp,i)
       Do j=1,i;       itr2=iptar(ilsp,j)
        ij=(i-1)*i/2+j; itr=(itr1-1)*itr1/2+itr2
        s = (tmar(ij)*tmar(ij)+tmai(ij)*tmai(ij))*g
        if(itr1.eq.itr2.and.i.ne.j) s = s + s
        om(itr) = om(itr) + s
       End do
      End do

      nom = iop(ie)*(iop(ie)+1)/2

      if(ktr.gt.0) then

      if(iop(ie).ne.np) Stop ' iop & np - ???'

      Do i=kopen+1,nopen;  itr1=iptar(ilsp,i)
       Do j=1,nj;          itr2=iptar(ilsp,j)
        ij=ntr+(i-kopen-1)*nj+j; itr=nom+(itr1-np-1)*ni+itr2
        s = (tmar(ij)*tmar(ij)+tmai(ij)*tmai(ij))*g
        if(itr1.eq.itr2.and.i.ne.j) s = s + s
        om(itr) = om(itr) + s
       End do
      End do

      end if

      nom = nom + ni * (jop(ie)-iop(ie))

      if(nom.gt.npl(ie)) then
       write(*,*) 'ilsp,nom,npl(ie)',ilsp,nom,npl(ie)
       write(*,'(a,10i5)') 'nopen,kopen,iop(ie),jop(ie),ktr', &
                            nopen,kopen,iop(ie),jop(ie),ktr

        Stop  'nom > npl(ie)'

      end if

      fl(1:nom,ie) = fl(1:nom,ie) + om(1:nom)

      write(nuq,'(F10.6,6i8,a)')  e(ie),ilsp,nom,iop(ie),jop(ie),np,ni, &
                              '   e(ie),ilsp,nom,iopen,jopen,np,ni'
      write(nuq,'(5D16.8)') (om(i),i=1,nom)

      go to 3
    4 Continue
      Close(nut)

!----------------------------------------------------------------------
! ... energy order:

      Allocate(ipe(ne));  Call SortR(ne,E,ipe)

! ... output omega:

      Do je=1,ne; ie=ipe(je)
       ntr = npl(ie)
       write(nuo,'(F10.6,5i8,a)')  e(ie),ntr,iop(ie),jop(ie),np,ni, &
                               '   e(ie),ntr,iopen,jopen,np,ni'
       write(nuo,'(5D16.8)') (fl(i,ie),i=1,ntr)
      End do
      Close (nuo)

! ... log information:

      Open(pri,file=AF_log)

      write(pri,*) 'Absent:'
      Do ie=1,ne; Do ilsp=1,nlsp
       if(ipl(ilsp,ie).gt.0) Cycle
       if(jopen(e(ie),ilsp).eq.0) Cycle
       write(pri,'(f10.6,i5)') E(ie),ilsp
      End do; End do

      write(pri,*) 'In excess:'
      Do ie=1,ne; Do ilsp=1,nlsp
       if(ipl(ilsp,ie).le.1) Cycle
       write(pri,'(f10.6,2i5)') E(ie),ilsp,ipl(ilsp,ie)
      End do; End do

      End  ! utility tmb_omb

