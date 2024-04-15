!======================================================================
!     UTILITY     sec_omn
!
!     target, zarm.om (zarm.omb, zarm.omb_top)  -> tr_###_###
!
! Call as   sec_omn  itr1=.. itr2=.. jtr1=.. jtr2=..  [label]
!
! lavel - additional information in the outout names
!----------------------------------------------------------------------

      Use target

      Implicit real(8) (a-h,o-z)

      Real(8), Allocatable ::  e(:), om(:), ecs(:), matr(:)
      Real(8), Allocatable ::  eom(:), ome(:)
      Integer, Allocatable ::  iom(:), ipt(:)

      Real(8), parameter :: Ry = 13.6056

      Character(80)  :: AS,AT,AF, label = ' ', sigma = 'sigma'
      Character(3) :: AA = 'tr_'

! ... files:

      Integer :: nut=1;   Character(20) :: targ   = 'target'
      Integer :: nuo=3;   Character(20) :: omega  = 'zarm.om'
                          Character(20) :: omegb  = 'zarm.omb'
                          Character(20) :: omegt  = 'zarm.omb_top'

      Integer :: out=9;   Character(80) :: AF_tr  = 'tr_nnn_nnn'

      Call inf_sec_omn
!----------------------------------------------------------------------
! ... target information:

      Call Check_file(targ)
      Open(nut,file=targ)
      Call R_target(nut)
      np=ntarg; Call Read_ipar(nut,'np',np)
      ni=ntarg; Call Read_ipar(nut,'ni',ni)
      Close(nut)

      ion = nz-nelc; zion=1.d0; if(ion.ne.0) zion=ion*ion
      E1 = etarg(1); etarg = (etarg-E1)*2.0
      Z = nz;  AWT = 0.d0

      Call Conv_au (Z,AWT,au_cm,au_eV,0)

! ... define units for output:

      i16 = 16;  Call Read_iarg('i16',i16)

! ... additional label if any:

      label = ' '
      Call Read_aarg('label',label)
      ilabel = LEN_TRIM(label)
      if(ilabel.gt.0) sigma = label

!-----------------------------------------------------------------------
! ... define transition under consideration:

      i11 = 1; Call Read_iarg('itr1',i11)
      i12 = 0; Call Read_iarg('itr2',i12);  if(i12.eq.0) i12=ntarg
      i21 = 1; Call Read_iarg('jtr1',i21)
      i22 = 0; Call Read_iarg('jtr2',i22);  if(i22.eq.0) i22=ntarg

      if(i12.gt.np) i12=np
!-----------------------------------------------------------------------
! ... read omega:

      if(Icheck_file(omegt).eq.1) then
       AF = omegt
      elseif(Icheck_file(omegb).eq.1) then
       AF = omegb
      elseif(Icheck_file(omega).eq.1) then
       AF = omega
      end if
      Call Read_aarg('om',AF)
      write(*,'(a,a)') 'Used  ',trim(AF)
      Call Check_file(AF)
      open(nut,file=AF)

      if(AF.eq.'zarm.omb_top') AA='cr_'


      mdim = ntarg*(ntarg+1)/2;  Allocate(matr(mdim))

      ip = 0; me = 0
    1 read(nut,*,end=2) ek,ns
      read(nut,*) (matr(i),i=1,ns)
      me = me +1
      ip = ip + ns
      go to 1
    2 Continue

      write(*,'(a,i10,a)') 'me = ',me,' - number of energies'
      write(*,'(a,i10,a)') 'ip = ',ip,' - number of omega entries'
      Allocate(ome(ip), eom(me), iom(0:me), ipt(me))
      write(*,'(a,f10.2,a)') 'memory = ', (ip*8 + me*8 +me*4)/1024.0/1024,  '  Mb'

      ip = 0; ne=0; iom=0
      rewind(nut)
    3 read(nut,*,end=4) ek,ns
      read(nut,*) (matr(i),i=1,ns)
      ne = ne + 1
      jp = ip + 1
      ip = ip + ns
      eom(ne) = ek
      iom(ne) = ip
      ome(jp:ip) = matr(1:ns)
      go to 3
    4 Close(nut)

      Call SortR(me,eom,IPT)

!-----------------------------------------------------------------------
! ... cycle over transitions:

      Allocate(e(0:me),om(0:me))

      Do itr1 = i11,i12
      Do itr2 = i21,i22

       if(itr1.gt.itr2) Cycle

       jtr1=itr1; jtr2=itr2

! ...  extract omega:

       itr = itr2*(itr2-1)/2 + itr1
       if(ion.ne.0) itr = itr - itr2 + 1

       om = 0.d0;  e(0) = etarg(itr2);  om(0)=0;  ne=0
       Do j=1,me; i=IPT(j)
        if(iom(i)-iom(i-1).lt.itr) Cycle
        ne = ne + 1
        e(ne)  = eom(i)
        om(ne) = ome(iom(i-1)+itr)
       End do

       de = Etarg(itr2)-Etarg(itr1)

! ... output:

      i1 = itr1; i2 = itr2

! ... define name for output file:

      write(AF_tr,'(a,i3.3,a,i3.3,a)') AA,i1,'_',i2
      if(ilabel.gt.0) &
      write(AF_tr,'(a,i3.3,a,i3.3,a,a,a)') AA,i1,'_',i2,'.',trim(label)

      Open(out,file=AF_tr)

! ... output data:

      de = etarg(i2)-etarg(i1)

! ... statistical weight for initial state:

      g=iabs(IStarg(i1))*(2.0*Ltarg(i1)+1)
      if(IStarg(i1).eq.0) g=Ltarg(i1)+1

      AT  =  'in ao^2'
      if(i16.gt.0) write(AT,'(a,i2.2)') 'in 10^-',i16

      if(ion.eq.0) then
       write(out,'(a,6x,a,10x,a,12x,a,10x,a,15x,a)') '#','eV',trim(sigma),'Ry','om',trim(AT)
      else
       write(out,'(a,6x,a,10x,a,12x,a,10x,a,14x,a,12x,a)') &
         '#','eV',trim(sigma),'Ry','om','Ry/z^2',trim(AT)
      end if

      if(ion.eq.0.and.i1.ne.i2) &
       write(out,'(f14.8,e16.8,f12.6,e16.8)') &
        etarg(i2)*Ry,0.d0,etarg(i2),0.d0

      Do ie=1,ne;   if(om(ie).eq.0.d0) Cycle
       es=e(ie)-etarg(i1)

       s=om(ie)/g/es                                 ! in [pi ao^2]
       if(i16.eq.0) s = s * 3.1415926                ! in ao^2
       if(i16.gt.0) s = s * 0.28003 * 3.1415926 * 10**(i16-16)   ! in 10^-16

       if(ion.eq.0) then
        write(out,'(f14.8,e16.8,f12.6,e16.8,f12.6)') &
                          e(ie)*Ry,s,e(ie),om(ie)
       else
        write(out,'(f14.8,e16.8,f12.6,e16.8,f12.6)') &
                          e(ie)*Ry,s,e(ie),om(ie),e(ie)/zion
       end if

      End do

      Close(out)

      End do; End do  ! over transitions

      End   !  sec_omn


!======================================================================
      Subroutine inf_sec_omn
!======================================================================
!     provide screen information about sec_omn utility
!----------------------------------------------------------------------

      Character :: A

      iarg = command_argument_count()
      if(iarg.eq.0) Return
      Call GET_COMMAND_ARGUMENT(1,A)
      if(A.ne.'?'.and.A.ne.'!') Return

      write(*,'(a)') &
'                                                                              ',&
'     "sec_omn" provides cross sections in table form                          ',&
'                                                                              ',&
'     zarm.om  + target  =>   tr_iii_fff[.label]                               ',&
'                                                                              ',&
'     additional possible input:  zarm.om, zarm_omb_par, om=..                 ',&
'                                                                              ',&
'     Arguments: itr1,itr2, jtr1,jtr2 - range for transition indexes           ',&
'                                       to be considered                       ',&
'                                       (itr -initial index, jtr - final)      ',&
'                                                                              ',&
'                i16 - defines units for output:                               ',&
'                      0 - sigma in a.u.  (default)                            ',&
'                      1 - sigma in 10-16 cm^2                                 ',&
'                      2 - sigma in 10-18 cm^2                                 ',&
'                     >2 - sigma in 10-i16 cm^2                                ',&
'                                                                              ',&
'                lavel - additional information in the outout names            ',&
'                                                                              ',&
'     Call as:   sec_omn itr1=.. itr2=.. jtr1=.. jtr2=.. [i16=.. label=..]     ',&
'                                                                              ',&
'     Results:   tr_###_###, where ### stand for transition indexes            ',&
'                                                                              '
      Stop ' '

      End Subroutine inf_sec_omn




