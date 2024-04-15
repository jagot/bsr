!=====================================================================
!     UTILITY  sec_om
!=====================================================================
!     provides collizion strengths or cross sections
!     for given transitions from information in 'zarm.om'
!
!     Arguments: i11,i12,  i21,i22  - range for first and second
!                                     transition indeces
!                i16 - control the output:
!                      0 - sigma in a.u.  (default)
!                      1 - sigma in 10-16
!                      2 - sigma in 10-18
!                     -1 - omega
!---------------------------------------------------------------------
      Use target

      Implicit real(8) (A-H,O-Z)

      Real, allocatable :: e(:), om(:), y(:), fl(:)
      Real :: Ry = 13.6057
      Logical :: EX
      Character(3)  :: AT
      Integer, allocatable :: itf1(:), itf2(:)

      Integer :: nut=1;  Character(20) :: AF_t  = 'target'
      Integer :: nuo=2;  Character(20) :: AF_om = 'zarm.om'
      Integer :: out=3;  Character(20) :: AF_tr = 'tr_##_##.dat'
      Integer :: nuf=8;  Character(20) :: AF_f  = 'f_values'
      Integer :: nur=9;  Character(20) :: AF_rr = 'rr_##_##.dat'

      Call inf_sec_om
!----------------------------------------------------------------------
! ... target information:

      Inquire(file=AF_t,EXIST=EX)
      if(.not.EX) then
       Stop  ' No target file: run H_targb ! '
      else
       Open(nut,file=AF_t); Call R_target(nut); close(nut)
      end if

      ion = nz-nelc; zion=1.d0; if(ion.gt.1) zion=ion*ion

      E1 = etarg(1); etarg = (etarg-E1)*2.0

      Z = nz;  AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,0)
      Ry = au_eV/2.d0

!----------------------------------------------------------------------
! ... f_values:

      ntf = 0
      if(Icheck_file(AF_f).eq.0) go to 25
      Open(nuf,file=AF_f)
      Call Read_ipar(nuf,'nt',ntf)
      if(ntf.eq.0) go to 25
      Allocate(itf1(ntf),itf2(ntf),fl(ntf))
      Do i = 1,ntf
       read(nuf,*) itf1(i),itf2(i),fl(i)
      End do
   25 Continue
      write(*,'(a,i6)') 'ntf =',ntf

!----------------------------------------------------------------------
! ... define transition under consideration:

      i11=1; i12=1; i21=1; i22=1
      iarg = IARGC()
      if(iarg.ge.1) then; Call GETARG(1,AT); read(AT,'(i3)') i11; end if
      if(iarg.ge.2) then; Call GETARG(2,AT); read(AT,'(i3)') i12; end if
      if(iarg.ge.3) then; Call GETARG(3,AT); read(AT,'(i3)') i21; end if
      if(iarg.ge.4) then; Call GETARG(4,AT); read(AT,'(i3)') i22; end if
      if(i12.lt.i11) i12=i11; if(i12.gt.ntarg) i12=ntarg
      if(i22.lt.i21) i22=i21; if(i22.gt.ntarg) i22=ntarg

      i16 = 0
      if(iarg.gt.4) then; Call GETARG(5,AT); read(AT,'(i3)') i16; end if

!----------------------------------------------------------------------
! ... define energies:

      Open(nuo,file=AF_om,status='OLD')

      me=0
      rewind(nuo)
    1 read(nuo,*,end=2) x,ns,(S,ix=1,ns)

      me = me + 1
      go to 1
    2 Continue
      write(*,*) ' me = ',me

! ... allocations:

      ntr = ntarg*(ntarg+1)/2;  Allocate(y(ntr),E(me),OM(me))

!----------------------------------------------------------------------
! ... cycle over transitions:

      Do i1=i11,i12; Do i2=i21,i22

       if(ion.ne.0.and.i1.eq.i2) Cycle

       if(i1.gt.i2) Cycle

! ... position of transition in arrays:

       itr=0; ij=0
       do i=1,ntarg
         jup=i; if(ion.gt.0) jup=i-1
         do j=1,jup
          ij=ij+1; if(i1.eq.j.and.i2.eq.i) itr=ij
         end do
       end do

! ... statistical weight for initial state:

       g=iabs(IStarg(i1))*(2.0*Ltarg(i1)+1)
       if(IStarg(i1).eq.0) g=Ltarg(i1)+1
       g1=iabs(IStarg(i1))*(2.0*Ltarg(i1)+1)
       if(IStarg(i1).eq.0) g1=Ltarg(i1)+1
       g2=iabs(IStarg(i2))*(2.0*Ltarg(i2)+1)
       if(IStarg(i2).eq.0) g2=Ltarg(i2)+1

! ... read data for given transition:

      ne=0
      rewind(nuo)
    5 read(nuo,*,end=15) x,ns,(y(i),i=1,ns)
      if(ns.lt.itr) go to 5
      if(x.le.etarg(i2)) go to 5

      ie=0
      Do i=1,ne; if(x.ne.e(i)) Cycle; ie=i; Exit; End do
      if(ie.eq.0) then
       ne=ne+1; if(ne.gt.me) Stop 'ne > me';  ie=ne
      end if
      e(ie)=x;  om(ie)=y(itr)

      go to 5

   15 if(ne.eq.0) Cycle

! ... ordering the data:

      Do i=1,ne-1; Do j=i+1,ne
       if(e(i).gt.e(j)) then
        S=e(i); e(i)=e(j); e(j)=S
        S=om(i); om(i)=om(j); om(j)=S
       end if
      End do; End do

! ... define name for output file:

      write(AF_tr,'(a,i3.3,a,i3.3)') 'tr_',i1,'_',i2

      Open(out,file=AF_tr)

! ... output data:

      de = etarg(i2)-etarg(i1)

      if(i16.eq.0) &
      write(out,'(a1,5a14,10x,a10)') '#','eV','sigma','Ry','om','Ry/z^2','in ao^2'
      if(i16.eq.1) &
      write(out,'(a1,5a14,10x,a10)') '#','eV','sigma','Ry','om','Ry/z^2','in 10^-16'
      if(i16.eq.2) &
      write(out,'(a1,5a14,10x,a10)') '#','eV','sigma','Ry','om','Ry/z^2','in 10^-18'

      if(i16.eq.-1) write(out,'(a1,5a14,10x,a10)') '#','Ry','om'
      if(i16.eq.-2) &
      write(out,'(a1,6a14,10x,a10)') '#','Ry','om','e1','s1','e2','s2','in 10^-16'


      if(ion.eq.0.and.i1.ne.i2.and.i16.ge.0) &
      write(out,'(f14.8,e16.8,f12.6,e16.8,f12.6)') etarg(i2)*Ry,0.d0,etarg(i2),0.d0,0.d0
      if(i16.eq.-1) write(out,'(f14.6,e16.8)') etarg(i2),0.d0

      if(i16.eq.-2) write(out,'(3(f14.6,e16.8))') &
         etarg(i2), 0.d0, etarg(i2)-etarg(i1), 0.d0, 0.d0

      Do ie=1,ne

       es=e(ie)-etarg(i1)
       s=om(ie)/g/es                                   ! in [pi ao^2]
       if(i16.eq.0) S = S * 3.1415926                  ! in ao^2
       if(i16.eq.1) s = s * 0.28003 * 3.1415926        ! in 10^-16
       if(i16.eq.2) s = s * 28.003  * 3.1415926        ! in 10^-18

       if(i16.ge.0) &
        write(out,'(f14.8,e16.8,f12.6,e16.8,f12.6)') &
                         e(ie)*Ry,s,e(ie),om(ie),e(ie)/zion
       if(i16.eq.-1) &
        write(out,'(f14.6,e16.8)') e(ie),om(ie)  ! only omega

       if(i16.eq.-2) then
        sn = 0.28003 * 3.1415926
        e1 = e(ie)-etarg(i1)
        e2 = e(ie)-etarg(i2)
        s1 = om(ie)/g1/e1 * sn
        s2 = om(ie)/g2/e2 * sn
        write(out,'(3(f14.6,e16.8))') e(ie),om(ie),e1,s1,e2,s2
       end if
      End do

      Close(out)

      if(i16.ne.-2) Cycle

      write(AF_rr,'(a,i3.3,a,i3.3)') 'rr_',i1,'_',i2
      Open(nur,file=AF_rr)

      f = 1.0
      Do i = 1, ntf
       if(itf1(i).ne.i1) Cycle
       if(itf2(i).ne.i2) Cycle
       f = fl(i)
       Exit
      End do

      write(nur,'(a,4f12.5)') 'f =',f,g,de

      f = 4 * g * f / de

      write(nur,'(a,f12.5)') 'norm =',f

      Do ie = 1,ne
       er = log( (e(ie)-etarg(i2))/de + 2.71828 )
       or = om(ie) / er / f
       write(nur,'(f14.8,3e16.8)') e(ie),om(ie),er,or
      End do

      End do; End do  ! over transitions

      End  ! UTILITY  sec_om



!======================================================================
      Subroutine inf_sec_om
!======================================================================
!     provide screen information about sec_om utility
!----------------------------------------------------------------------

      Implicit real(8) (a-h,o-z)

      Character(80) :: A

      Call get_command_argument(1,A)
      if(A.ne.'?') Return

      write(*,*) &
'                                                                ',&
'     SEC_OM provides cross sections in table form               ',&
'                                                                ',&
'     zarm.om + target  =>   tr_ii_ff.dat                        ',&
'                                                                ',&
'     Call as:  sec_om   ii1 ii2 ff1 ff2 [i16]                   ',&
'                                                                ',&
'     ii1,ii2 - range for initial index ii                       ',&
'     ff1,ff2 - range for final   index ff                       ',&
'                                                                ',&
'     i16 - control the output:                                  ',&
'           0 - sigma in a_o^2 (default)                         ',&
'           1 - sigma in 10-16 cm^2                              ',&
'           2 - sigma in 10-18 cm^2                              ',&
'          -1 - only omega                                       ',&
'                                                                '
      Stop ' '

      End Subroutine inf_sec_om

