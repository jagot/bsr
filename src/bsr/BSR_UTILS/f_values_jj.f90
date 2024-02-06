!=======================================================================
!     Utility f-values  (h.nnn, target_jj -->  f_values_nnn
!=======================================================================
      Use target_jj
      Use channel_jj

      Implicit real(8) (A-H,O-Z)
      Real(8), allocatable :: CF(:,:,:)
      Character(80) :: AF

      Call get_command_argument(1,AF)
      if(AF.eq.'?') then    !  help section
       write(*,*)
       write(*,*) 'f_values_jj calculates the f-values on the base of asymptotic'
       write(*,*)
       write(*,*) 'coefficients stored in file h.nnn'
       write(*,*)
       write(*,*) 'h.nnn, target_jj  -->  f_values_nnn'
       write(*,*)
       write(*,*) 'Call as:   f_values [klsp=..]'
       write(*,*)
       write(*,*) 'klsp [1] - the partial-wave index, nnn'
       write(*,*)
       write(*,*) 'if h.nnn is absent, utility will check H.DAT file'
       write(*,*)
       Stop ' '
      end if

      klsp=1; Call Read_iarg('klsp',klsp)

! ... target information:

      nut = 1
      Open(nut,file='target_jj',status='OLD')
      Call Read_target_jj(nut)
      Call Read_channel_jj(nut,klsp)
      close(nut)

! ... read  assymtotic coefficients from H.DAT file

      in = 2
      write(AF,'(a,i3.3)') 'h.',klsp
      if(Icheck_file(AF).eq.0) AF = 'H.DAT'
      Call Read_aarg('h',AF)
      Open(in,file=AF,status='OLD',form='UNFORMATTED')

      read(in) NELC, NZ, LRANG2, LAMAX, NAST, RA, BSTO

      read(in) (E,N=1,NAST)
      read(in) (L,N=1,NAST)
      read(in) (I,N=1,NAST)
      read(in) ((C,K=1,3),L=1,LRANG2)

      read(in) LRGL, NSPN, NPTY, NCHAN, MNP2, MORE0

      read(in) (NCONAT, N=1,NAST)
      read(in) (L2P, I=1,NCHAN)

      Allocate (CF(nch,nch,lamax))

      read(in) (((CF(I,J,K), I=1,NCH), J=1,NCH), K=1,LAMAX)

! ... calculations of f-values:

      nu = 3
      write(AF,'(a,i3.3)') 'f_values_',klsp
      Open(nu,file=AF)

      eps = 1.D-6
      Call f_values(nu,eps,lamax,CF)

      End ! program


!======================================================================
      Subroutine f_values(pri,eps_acf,km,ACF)
!======================================================================
!     define the f-values between target states based on the
!     given asimptotic coefficients ACF for k=1
!----------------------------------------------------------------------
      Use zconst, only: c_au, time_au
      Use target_jj
      Use channel_jj

      Implicit none
      Real(8) :: ACF(nch,nch,km)
      Integer :: pri,km
      Real(8) :: eps_acf

      Real(8), allocatable :: AK(:,:), AS(:,:)
      Integer, allocatable :: IP(:,:)
      Real(8) :: S,SS, g1,g2, de, a,f
      Integer :: i,j, i1,i2,nt
      Real(8), external :: CLEBCH, Z_6jj, Reduce_factor

      if(.not.allocated(AK)) Allocate(AK(ntarg,ntarg)); AK=0.d0
      if(.not.allocated(AS)) Allocate(AS(ntarg,ntarg)); AS=0.d0
      if(.not.allocated(IP)) Allocate(IP(ntarg,ntarg)); IP=0

      Do i=1,nch-1; Do j=i+1,nch

       i1=iptar(i); i2=iptar(j);   if(i1.eq.i2) Cycle
       S = ACF(i,j,1)/2.d0;        if(abs(S).lt.eps_acf) Cycle
       SS = Reduce_factor(i,j,1);  if(abs(SS).lt.eps_acf) Cycle

       de=Etarg(i2)-Etarg(i1)

       S = S / SS

       S = S*S
       g1 = jtarg(i1)+1
       g2 = jtarg(i2)+1

       f = 2.d0/3.d0*de*S /g1
       a = 4.d0/3.d0*de**3*S/c_au**3/time_au /g2

       AS(i1,i2) = AS(i1,i2) + S
       AK(i1,i2) = AK(i1,i2) + f
       AK(i2,i1) = AK(i2,i1) + a
       IP(i1,i2) = IP(i1,i2) + 1

      End do; End do

      nt = 0
      Do i1=1,ntarg; Do i2=i1,ntarg
       if(IP(i1,i2).eq.0) Cycle
       AS(i1,i2) = AS(i1,i2) / IP(i1,i2)
       AK(i1,i2) = AK(i1,i2) / IP(i1,i2)
       AK(i2,i1) = AK(i2,i1) / IP(i1,i2)
       nt = nt + 1
      End do; End do

! ... total decay probabilities:

      Do i=2,ntarg;  AK(i,i) = SUM(AK(i,1:i-1)); End do

! ... print results:

      if(nt.gt.0) &
      write(pri,'(i10,3x,a,a,i8)') nt,'Target radiative data:  ', &
        '  s-value, f-value, A-value, branching ratio, nt =',nt

      Do i=1,ntarg-1; Do j=i+1,ntarg
       if(AK(i,j).eq.0.d0) Cycle

       write(pri,'(2i5,1PE12.3,2E12.3,0Pf10.5,5x,a,5x,a)') &
                   i,j,AS(i,j),AK(i,j),AK(j,i),AK(j,i)/AK(j,j), &
                   trim(AFT(i)),trim(AFT(j))
      End do; End do

! ... polarizability of the ground state:

      a = 0.d0
      Do i = 2,ntarg
       de = Etarg(i)-Etarg(1)
       a = a + AK(1,i)/(de*de)
      End do

      if(a.ne.0.d0) &
      write(pri,'(/a,f10.3/)') 'Polarizability of the ground state =',a

      End Subroutine f_values


!======================================================================
      Real(8) Function Reduce_factor(ich,jch,k)
!======================================================================
!     define factor connecing reduced dipole matrix element with
!     asymptotic coefficient for multipole index k
!----------------------------------------------------------------------
      Use target_jj
      Use channel_jj

      Implicit none
      Integer, intent(in) :: ich,jch,k
      Integer :: j1,j2,JT1,JT2,JJ,kz
      Real(8) :: S, zero = 0.d0
      Real(8), external :: Cjkj, Z_6j2
      Integer, external :: j_kappa

      Reduce_factor = zero

      j1 = j_kappa(kch(ich))
      j2 = j_kappa(kch(jch))
      S = Cjkj(j1,k,j2)
      if(S.eq.zero) Return
      JT1 = jtarg(iptar(ich))
      JT2 = jtarg(iptar(jch))
      JJ  = jpar

      S = S * Z_6j2(JT1,j1,JJ,j2,JT2,k+k)
      kz = JT2+j1+JJ; kz = kz/2
      S = S * (-1)**kz

      Reduce_factor = S

      End Function Reduce_factor
