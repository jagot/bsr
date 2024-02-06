!----------------------------------------------------------------------
!     h.001 + h.002 + ... + h.nnn  --->  H.DAT
!----------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)

      Real(8),allocatable :: ENAT(:), ENAT1(:), VALUE(:)
      Integer,allocatable :: LAT(:), LAT1(:), ISAT(:), ISAT1(:), &
                             NCONAT(:), L2P(:), IPAT(:), IPAT1(:), jkch(:)
      Real(8),allocatable :: COEFF(:,:), WMAT(:,:), CF(:,:,:)

      Integer :: ih1    = 1
      Integer :: ih2    = 999
      Integer :: MRANG2 = 201
      Integer :: folder = 0

      Character(80) :: AF

! ... Check the arguments:

      Call get_command_argument(1,AF)
      if(AF.eq.'?') then    !  help section
       write(*,*)
       write(*,*) 'sum_hh merges the different h.nnn files: '
       write(*,*)
       write(*,*) 'h.001 + h.002 + ... + h.nnn --->  H.DAT'
       write(*,*)
       write(*,*) 'Call as:   sum_hh  [ih1=.. ih2=.. folder=..]'
       write(*,*)
       write(*,*) 'with default values for partial h.nnn: ih1=1, ih2=999'
       write(*,*)
       write(*,*) 'if folder <> 0, we supposed paths:  nnn/h.nnn'
       write(*,*)
       Stop ' '
      end if

! ... files:

      in=1
      iout=2; Open(iout,file='H.DAT',form='UNFORMATTED')
      ipri=9; Open(ipri,file='sum_hh.log')


      Call Read_iarg('ih1',ih1)
      Call Read_iarg('ih2',ih2)

! ... Cycle over different H-files:

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

        Allocate (COEFF(3,MRANG2), ENAT(NAST),LAT(NAST),ISAT(NAST), &
                  ENAT1(NAST), LAT1(NAST), ISAT1(NAST),NCONAT(NAST), &
                  IPAT(NAST), IPAT1(NAST) )

        read(in) ENAT
        read(in) LAT
        read(in) ISAT,IPAT

        COEFF = 0.d0
        LM = min(MRANG2,LRANG2)
        read(in) ((COEFF(K,L),K=1,3),L=1,LM)

        NELC1=NELC; NZ1=NZ; NAST1=NAST; ENAT1=ENAT; LAT1=LAT; ISAT1=ISAT; IPAT1=IPAT
        LAMAX1=LAMAX; RA1=RA; BSTO1=BSTO

        write(iout) NELC, NZ, MRANG2, LAMAX, NAST, RA, BSTO
        write(iout) ENAT
        write(iout) LAT
        write(iout) ISAT,IPAT
        write(iout) ((COEFF(K,L),K=1,3),L=1,MRANG2)
        istart = 0

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

        read(in) ENAT
        read(in) LAT
        read(in) ISAT,IPAT
        LM = min(MRANG2,LRANG2)
        read(in) ((COEFF(K,L),K=1,3),L=1,LM)

        Do N=1,NAST
         if(ENAT1(N).ne.ENAT(N)) then
          write(*,*) ' ENAT <> ENAT1 for ',AF; Stop
         end if
         if(LAT1(N).ne.LAT(N)) then
          write(*,*) ' LAT <> LAT1 for ',AF; Stop
         end if
         if(ISAT1(N).ne.ISAT(N)) then
          write(*,*) ' ISAT <> ISAT1 for ',AF; Stop
         end if
         if(IPAT1(N).ne.IPAT(N)) then
          write(*,*) ' IPAT <> IPAT1 for ',AF; Stop
         end if
        End do

      end if

!----------------------------------------------------------------------
!                                               scattering information:

    1 read(in) LRGL, NSPN, NPTY, NCHAN, MNP2, MORE0

      MORE=1; if(ih.eq.ih2.and.MORE0.eq.0) MORE=0

      write(ipri,'(a,3i4,a,i5,a,i6)') 'LSP=', LRGL, NSPN, NPTY, &
        '   nchan =',nchan,'   nhm =',MNP2

      Allocate (L2P(nchan), CF(nchan,nchan,lamax), VALUE(MNP2),&
                WMAT(NCHAN,MNP2),jkch(nchan))

      read(in) (NCONAT(N), N=1,NAST)
      read(in) (L2P(I), I=1,NCHAN),(jkch(i),i=1,nchan)
      read(in) (((CF(I,J,K), I=1,NCHAN), J=1,NCHAN), K=1,LAMAX)
      read(in) (VALUE(K),K=1,MNP2)
      read(in) ((WMAT(I,K),I=1,NCHAN), K=1,MNP2)

      write(iout) LRGL, NSPN, NPTY, NCHAN, MNP2, MORE
      write(iout) (NCONAT(N), N=1,NAST)
      write(iout) (L2P(I), I=1,NCHAN),(jkch(i),i=1,nchan)
      write(iout) (((CF(I,J,K), I=1,NCHAN), J=1,NCHAN), K=1,LAMAX)
      write(iout) (VALUE(K),K=1,MNP2)
      write(iout) ((WMAT(I,K),I=1,NCHAN), K=1,MNP2)

      write(ipri,'(5f14.5)') VALUE(mnp2-4:mnp2)

      Deallocate(L2P, CF, VALUE, WMAT, jkch)

      istart = 0
      if(MORE0.ne.0) go to 1

      End do  ! over ih

      End  ! sum_hh
