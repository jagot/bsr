!======================================================================
!     H.DAT  -->  target  (bsr-format)
!======================================================================

      Use target, nconat => ictarg
      Use channels

      Implicit real(8) (A-H,O-Z)

      Allocatable  llch(:), jptar(:), eval(:), aval(:,:)

      Character(10) :: AF = 'H.DAT'
      Character(3) :: AT = '000'
      Character(1), external :: AL

      Real(8), parameter ::  Ry = 13.6057

!----------------------------------------------------------------------
! ... files:

      iarg = IARGC(); if(iarg.gt.0)  Call GETARG(1,AF)
      if(AF.eq.'?') then
        write(*,'(/a)') 'h_targb creates "target" file based on the "H.DAT" file'
        write(*,'(/a)') '     H.DAT  -->  target  (bsr-format)'
        write(*,'(/a)') 'Call as:  h_targb [name]'
        write(*,'(/a)') 'with default name - H.DAT'
        Stop ' '
      end if

      Call Check_file(AF)
      in=1;   Open(in,file=AF,form='UNFORMATTED',status='OLD')
      isch=2; Open(isch,status='SCRATCH',form='UNFORMATTED')
      iout=3; Open(iout,file='target')

!----------------------------------------------------------------------
! ... target states information:

      read(in) NELC, NZ, LRANG2, LAMAX, ntarg, RA, BSTO

      m=ntarg; Call Allocate_target(m)

      read(in) etarg
      read(in) ltarg

      iptarg = 1
      read(in,err=10) istarg,iptarg
      go to 20

   10 rewind(in)
      read(in) NELC, NZ, LRANG2, LAMAX, ntarg, RA, BSTO
      read(in) etarg
      read(in) ltarg
      read(in) istarg
   20 Continue

! ... Buttle corrections:

      read(in) ((C,k=1,3),l=1,LRANG2)

      coupling = 'LS'; if(istarg(1).eq.0) coupling= 'JK'

!----------------------------------------------------------------------
! ... read symmetries and record useful information in scratch file:

      nlsp=0; mch=0

    1 read(in,end=2) LRGL, NSPN, NPTY, NCHAN, MNP2, MORE

      nlsp=nlsp+1; if(nchan.gt.mch) mch=nchan

      read(in) NCONAT

      Allocate(llch(NCHAN), jptar(NCHAN), eval(MNP2))

      read(in) llch

! ... asymptotic coefficients:

      read(in) (((C,i=1,NCHAN),j=1,NCHAN),K=1,LAMAX)

! ... (N+1)-electron eigenvalues:

      read(in) eval

! ... surface amplitudes:

      read(in) ((C,i=1,NCHAN),j=1,MNP2)

! ... find pointer 'channel --> target'

      i1=1
      Do it=1,ntarg
       if(NCONAT(it).eq.0) Cycle
       i2=i1+NCONAT(it)-1
       Do i=i1,i2; jptar(i)=it; End do
       i1=i2+1
      End do

! ... record for a while:

      write(isch) LRGL, NSPN, NPTY, NCHAN
      write(isch) llch
      write(isch) jptar
      write(isch) (eval(i),i=MNP2,MNP2-4,-1)

      Deallocate (llch,jptar,eval)

      if(more.eq.1) go to 1

!----------------------------------------------------------------------
! ... read the whole information from scratch file:

   2  Call Allocate_channels
      Allocate(aval(nlsp,5))
      rewind(isch)
      Do i=1,nlsp
       read(isch) lpar(i), ispar(i), ipar(i), nch(i)
       read(isch) (lch(i,j),j=1,nch(i))
       read(isch) (iptar(i,j),j=1,nch(i))
       read(isch) (aval(i,j),j=1,5)
      End do

!----------------------------------------------------------------------
! ... find jk-values in JK-case:

      jkch=0
      if(coupling.eq.'JK') then
       Do i = 1,nlsp
        K_min=iabs(lpar(i)-1)+1
        K_max=iabs(lpar(i)+1)+1
        Do j = 1,nch(i)
         ll = 2*lch(i,j)+1; JT = ltarg(iptar(i,j))+1
         Do KK = K_min,K_max,2
          if(ITRI(JT,ll,kk).eq.0) Cycle
          m = 0
          Do jj = 1,j-1
           if(lch(i,j).eq.lch(i,jj).and. &
              iptar(i,j).eq.iptar(i,jj).and. &
              jkch(i,jj).eq.kk) m=1
           if(m.eq.1) Exit
          End do
          if(m.eq.1) Cycle
          jkch(i,j) = kk; Exit
         End do
        End do
       End do
      end if

!----------------------------------------------------------------------
! ... record 'target' file:

      write(iout,'(a)') '   e + ...   '
      write(iout,'(72(''-''))')

      write(iout,'(a,a2,a,a2,a)') 'coupling = ',coupling,'    !   ', &
                                   coupling,' - case'
      write(iout,'(a,i3,a)') 'nz = ',nz,'         !   nuclear charge'
      write(iout,'(a,i3,a)') 'nelc =',nelc, &
                             '        !   number of electrons'
      write(iout,'(a,f8.3,a)') 'RA = ',RA,'    !   bouder radius'
      write(iout,'(a,i2,a)') 'lamax = ',lamax, &
                             '       !   max. multipole index'
      write(iout,'(72(''-''))')

      write(iout,'(a,i3,a)') 'ntarg =',ntarg, &
                             '       !   number of target states'
      write(iout,'(72(''-''))')
      Do i = 1,ntarg
       write(AF,'(a,i3.3)') 'targ_',i
       E_Ry = 2*(Etarg(i)-Etarg(1))
       E_eV = E_Ry * Ry
       write(iout,'(a1,2x,a10,3i5,F16.8,2i4,F10.6,f10.3)') &
         '?',AF,ltarg(i),istarg(i),iptarg(i),Etarg(i),0,0,E_Ry,E_eV
      End do
      write(iout,'(72(''-''))')

      write(iout,'(a,i3,5x,a)') 'nlsp =',nlsp, &
                             '   !   number of partial waves'
      write(iout,'(72(''-''))')
      Do i = 1,nlsp
       AF = 'no'
       i3=1-2*ipar(i)
       write(iout,'(i3.3,3i5,5x,2a10,2i5)') &
                    i,lpar(i),ispar(i),i3,AF,AF,0,0
      End do
      write(iout,'(72(''-''))')

      write(iout,'(a)') 'channels:'
      write(iout,'(72(''-''))')

      Do i = 1,nlsp
       write(iout,'(i3,a,i3.3,a,i4,a,2i6)')  &
                    i,'.  ',i,'  nch =',nch(i),'  nc =',0,0
       Do j = 1,nch(i)
        write(iout,'(2x,a1,a1,3i5,i10,i5)') &
              'k',AL(lch(i,j),1),lch(i,j),iptar(i,j),j,0,jkch(i,j)
       End do
       write(iout,'(72(''-''))')
      End do

      write(iout,'(a,i5)') 'max_ch =',mch
      write(iout,'(72(''-''))')

      End  ! utility H_targb
