!----------------------------------------------------------------------
!     h.001 + h.002 + ... + h.nnn --->  H.DAT    (JJ-coupling)
!----------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)

      Real(8),allocatable :: ENAT(:), ENAT1(:), VALUE(:)
      Integer,allocatable :: LAT(:), LAT1(:), ISAT(:), ISAT1(:), &
                             NCONAT(:), L2P(:), K2P(:), PAT(:), jptar(:)
      Real(8),allocatable :: COEFF(:,:), WMAT(:,:), CF(:,:,:)

      Integer, allocatable :: jpar(:), ipar(:),  nch(:), lch(:,:), kch(:,:), iptar(:,:)

      Integer :: ih1 = 1
      Integer :: ih2 = 999
      Integer :: MRANG2 = 201
      Integer :: folder = 0

      Character(80) :: AF, AS
      Character(1), external :: AL

! ... files:

      in=1
      iout=2; Open(iout,file='H.DAT',form='UNFORMATTED')
      isch=3; Open(isch,status='SCRATCH',form='UNFORMATTED')
      ipri=9; Open(ipri,file='sum_hh.log')

! ... Check the arguments:

      Call get_command_argument(1,AF)
      if(AF.eq.'?') then    !  help section
       write(*,*)
       write(*,*) 'sum_hh_jj merges the different h.nnn files after DBSR calculations: '
       write(*,*)
       write(*,*) 'Call as:   sum_hh  [ih1=.. ih2=.. folder=..]'
       write(*,*)
       write(*,*) 'with default values for partial h.nnn: ih1=1, ih2=999'
       write(*,*)
       write(*,*) 'folder[0] <> 0 supposed paths:  nnn/h.nnn'
       write(*,*)
       Stop ' '
      end if

      Call Read_iarg('ih1',ih1)
      Call Read_iarg('ih2',ih2)

! ... Cycle over different H-files:

      nlsp = 0; mch=0
      istart= 1
      Do ih = ih1,ih2

       write(AF,'(a,i3.3)') 'h.',ih
       if(folder.ne.0) write(AF,'(i3.3,a,i3.3)') ih,'/h.',ih

       if(Icheck_file(AF).eq.0) Cycle
       write(*,*) trim(AF)
       write(ipri,'(a)') trim(AF)
       Open(in,file=AF,form='UNFORMATTED',status='OLD')

!----------------------------------------------------------------------
!                                                   target information:

       read(in) NELC, NZ, LRANG2, LAMAX, NAST, RA, BSTO

       if(istart.eq.1)  then

        Allocate (COEFF(3,MRANG2), ENAT(NAST),LAT(NAST),ISAT(NAST),PAT(NAST), &
                  ENAT1(NAST), LAT1(NAST), ISAT1(NAST),  NCONAT(NAST))
        read(in) ENAT
        read(in) LAT             ! LAT --> JJ-value
        read(in) ISAT            ! ISAT --> 0

        COEFF = 0.d0
        LM = min(MRANG2,LRANG2)
        read(in) ((COEFF(K,L),K=1,3),L=1,LM)

        NELC1=NELC; NZ1=NZ; NAST1=NAST; ENAT1=ENAT; LAT1=LAT; ISAT1=ISAT
        LAMAX1=LAMAX; RA1=RA; BSTO1=BSTO

        write(iout) NELC, NZ, MRANG2, LAMAX, NAST, RA, BSTO
        write(iout) ENAT
        write(iout) LAT
        write(iout) ISAT
        write(iout) ((COEFF(K,L),K=1,3),L=1,MRANG2)
        istart=0

       else

        if(NAST1.ne.NAST) then
          write(*,*) ' NAST <> NAST1 for  ',AF; Stop
        end if
        if(NZ1.ne.NZ) then
          write(*,*) ' NZ <> NZ1 for  ',AF; Stop
        end if
        if(LAMAX1.ne.LAMAX) then
          write(*,*) ' LAMAX <> LAMAX1 for  ',AF; Stop
        end if
        if(RA1.ne.RA) then
          write(*,*) ' RA <> RA1 for  ',AF; Stop
        end if
        if(BSTO1.ne.BSTO) then
          write(*,*) ' BSTO <> BSTO1 for  ',AF; Stop
        end if

        read(in) ENAT1
        read(in) LAT1
        read(in) ISAT1
        LM = min(MRANG2,LRANG2)
        read(in) ((COEFF(K,L),K=1,3),L=1,LM)

        Do N=1,NAST
         if(ENAT1(N).ne.ENAT(N)) then
          write(*,*) ' ENAT <> ENAT1 for  ',AF; Stop
         end if
         if(LAT1(N).ne.LAT(N)) then
          write(*,*) ' LAT <> LAT1 for  ',AF; Stop
         end if
        End do

      end if

!----------------------------------------------------------------------
!                                               scattering information:

    1 read(in) LRGL, NSPN, NPTY, NCHAN, MNP2, MORE0

      MORE=1; if(ih.eq.ih2.and.MORE0.eq.0) MORE=0

      write(ipri,'(a,3i4,a,i5,a,i6)') 'LSP=', LRGL, NSPN, NPTY, &
        '   nchan =',nchan,'   nhm =',MNP2

      Allocate (L2P(nchan), K2P(nchan), jptar(nchan), CF(nchan,nchan,lamax), VALUE(MNP2),&
                WMAT(NCHAN,MNP2))

      read(in) (NCONAT(N), N=1,NAST)
      read(in) L2P, K2P
      read(in) CF
      read(in) VALUE
      read(in) WMAT

      write(iout) LRGL, NSPN, NPTY, NCHAN, MNP2, MORE
      write(iout) (NCONAT(N), N=1,NAST)
      write(iout) L2P, K2P
      write(iout) CF
      write(iout) VALUE
      write(iout) WMAT

      write(ipri,'(5f14.5)') VALUE(mnp2-4:mnp2)

      jptar = 0
      i1=1
      Do it=1,NAST
       if(NCONAT(it).eq.0) Cycle
       i2=i1+NCONAT(it)-1
       Do i=i1,i2; jptar(i)=it; End do
       i1=i2+1
      End do

      nlsp=nlsp+1; if(nchan.gt.mch) mch=nchan

      write(isch) LRGL, NPTY, NCHAN
      write(isch) L2P
      write(isch) K2P
      write(isch) jptar

      Deallocate(L2P, K2P, jptar, CF, VALUE, WMAT)

      if(MORE0.ne.0) go to 1

      End do

      Close(in)  ! over ih

!----------------------------------------------------------------------
! ... re-load the whole data and determine the parity of target states:

      Allocate(jpar(nlsp), ipar(nlsp), nch(nlsp),  &
               lch(nlsp,mch), kch(nlsp,mch), iptar(nlsp,mch))

      rewind(isch)
      Do i=1,nlsp
       read(isch) jpar(i), ipar(i), nch(i)
       read(isch) (lch(i,j),j=1,nch(i))
       read(isch) (kch(i,j),j=1,nch(i))
       read(isch) (iptar(i,j),j=1,nch(i))
      End do

      PAT = 0
      Do i=1,nlsp
       Do ich=1,nch(i); it=iptar(i,ich)
        ip = ipar(i);
        if(mod(lch(i,ich),2).eq.1) then
         if(ip.eq.0) then; ip=1; else; ip=0; end if
        end if
        PAT(it) = 1-2*ip
       End do
      End do

!----------------------------------------------------------------------
! ... record 'target' file:

      ires=10; Open(ires,file='target.hhh')

      ntarg = NAST

      write(ires,'(a)') '   e + ...   '
      write(ires,'(80(''-''))')

      write(ires,'(a,a2)') 'coupling = ','JJ'

      write(ires,'(a,i3,T19,a)') 'nz    =',nz,    '!   nuclear charge'

      write(ires,'(a,i3,T19,a)') 'nelc  =',nelc,  '!   number of electrons'
      write(ires,'(80(''-''))')

      write(ires,'(a,i3,T19,a)') 'ntarg =',ntarg, '!   number of target states'
      write(ires,'(80(''-''))')

      Do i = 1,ntarg

       write(AS,'(a,i3.3)') 'targ_',i
       write(AF,'(a,i3.3)') 't',i

       E_Ry = 2 * (ENAT(i)-ENAT(1))
       E_ev = E_Ry * 27.211269 / 2

       write(ires,'(a,7x,a,3i5,F16.8,2i4,F12.6,F10.5)') &
         trim(AF),trim(AS), LAT(i), 0, PAT(i), ENAT(i), 0,0, E_Ry,E_eV

      End do
      write(ires,'(80(''-''))')

      write(ires,'(a,i3,5x,a)') 'nlsp =',nlsp, &
                             '   !   number of partial waves'
      write(ires,'(80(''-''))')

      Do i = 1,nlsp
       ii=1-2*ipar(i)
       write(ires,'(i3.3,3i5,5x,2a10,2i5)')  i,jpar(i),0,ii   ! ,'no','no',0,0
      End do
      write(ires,'(80(''-''))')

      write(ires,'(a)') 'channels:'
      write(ires,'(80(''-''))')

      Do ilsp = 1,nlsp
       write(ires,'(i3,a,i3.3,a,i4,a,2i6)')  ilsp,'.  ',ilsp,'  nch =',nch(ilsp), &
            '  nc =',0,0
       Do i = 1,nch(ilsp)
        write(ires,'(2x,a1,a1,3i5,i10,i5)') &
              'k',AL(lch(ilsp,i),1),lch(ilsp,i),iptar(ilsp,i),i,0,j_kappa(kch(ilsp,i))
       End do
       write(ires,'(80(''-''))')
      End do

      write(ires,'(a,i5,T20,a)')   'max_ch =',mch,   '!  max. number of channels'
      write(ires,'(a,f8.3,T20,a)') 'RA = ',RA,       '!  bouder radius'
      write(ires,'(a,i2,T20,a)')   'lamax = ',lamax, '!  max. multipole index'


!----------------------------------------------------------------------
! ... record 'target_jj' file:

      nut=12; Open(nut,file='target_jj.hhh')

      write(nut,'(a)') '  e + ...'
      write(nut,'(80(''-''))')
      write(nut,'(a,i4,5x,a)') &
                'nz    = ',nz,   ' !   nuclear charge'
      write(nut,'(a,i4,5x,a)') &
                'nelc  = ',nelc, ' !   number of electrons'
      write(nut,'(80(''-''))')
      write(nut,'(a,i4,5x,a)') &
                'ntarg = ',ntarg,' !   number of target states'
      write(nut,'(80(''-''))')

      Do i=1,ntarg

       write(AS,'(a,i3.3)') 'targ_',i
       write(AF,'(a,i3.3)') 't',i

       E_Ry = 2 * (ENAT(i)-ENAT(1))
       E_ev = E_Ry * 27.211269 / 2

       write(nut,'(a,10x,a,2x,2i4,F18.8,2i4,F12.6,F10.5)') &
         trim(AF),trim(AS), LAT(i), PAT(i), ENAT(i), 0,0, E_Ry,E_eV

      End do

      write(nut,'(80(''-''))')

! ... partial waves:

      write(nut,'(a,i4,5x,a)') 'nlsp  = ',nlsp,' !   number of partial waves'
      write(nut,'(80(''-''))')

      Do i = 1,nlsp
       write(nut,'(i3,a,2i4)') i,'.',jpar(i),1-2*ipar(i)  ! ,'no','no',0,0
      End do
      write(nut,'(80(''-''))')

! ... channels:

      write(nut,'(a,i5)') 'channels:'
      write(nut,'(80(''-''))')

      Do ilsp = 1,nlsp

       write(nut,'(i3,a,i3.3,a,i4,a,2i6)')  ilsp,'.  ',ilsp,'  nch =',nch(ilsp), &
            '  nc =',0,0
       write(nut,'(80(''-''))')

       Do i = 1,nch(ilsp)
        write(nut,'(i3,a1,2x,a1,a1,2x,2i6,i8)') &
         i,'.','k',AL(lch(ilsp,i),1),kch(ilsp,i),iptar(ilsp,i),0
       End do
       write(nut,'(80(''-''))')

      End do ! ilsp

      write(nut,'(a,i5,T20,a)')   'max_ch =',mch,   '!  max. number of channels'
      write(nut,'(a,f8.3,T20,a)') 'RA = ',RA,       '!  bouder radius'
      write(nut,'(a,i2,T20,a)')   'lamax = ',lamax, '!  max. multipole index'
      write(nut,'(80(''-''))')

      End  ! sum_hh_jj
