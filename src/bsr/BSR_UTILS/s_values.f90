!=======================================================================
!     Program s-values
!=======================================================================
!     generate the oscillator strengths (dipole only)
!     based on the asymptotic coefficients saved in the H.DAT file
!
!     target,  H.DAT   ->   s_values,  s_values_###
!-----------------------------------------------------------------------
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)
      Character(80) :: AF
      Real(8), allocatable :: CF(:,:,:), OS(:,:)

!----------------------------------------------------------------------
      iarg = COMMAND_ARGUMENT_COUNT()
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)

      if(AF.eq.'?') then
        write(*,*)
        write(*,'(a)') 'generation dipole mtrix elements for transition between target states'
        write(*,'(a)') 'based on the asymptotic coefficients saved in the H.DAT file'
        write(*,*)
        write(*,'(a)') 'target,  H.DAT   ->   s_values,  s_values_###T'
        write(*,*)
        write(*,'(a)') 'Call as: s_values klsp1=... klsp2=... L1=... L2=...  h=...'
        write(*,*)
        write(*,'(a)') 'klsp1, klsp2 - restrictons on partial waves'
        write(*,'(a)') 'L1, L2 - restriction on total orbital moments'
        write(*,'(a)') 'h  -  alternative name for H.DAT file'
        Stop
      end if
!----------------------------------------------------------------------
! ... partial waves restrictions:

      klsp1=0;  Call Read_iarg('klsp1',klsp1)
      klsp2=0;  Call Read_iarg('klsp2',klsp2)

      LPAR1 = -1;  Call Read_iarg('L1',LPAR1)
      LPAR2 = -1;  Call Read_iarg('L2',LPAR2)

! ... target information:

      nut = 1
      Open(nut,file='target',status='OLD')
      Call R_target(nut)
      Call R_channels(nut)

      np=ntarg; Call Read_ipar(nut,'np',np)
      ni=ntarg; Call Read_ipar(nut,'ni',ni)

      Allocate(OS(ntarg,ntarg)); OS = 0.d0

! ... read  assymtotic coefficients from H.DAT file

      in = 2
      AF = 'h.dat'
      if(Icheck_file(AF).eq.0) AF = 'H.DAT'
      Call Read_aarg('h',AF)
      Open(in,file=AF,status='OLD',form='UNFORMATTED')

      read(in) NELC, NZ, LRANG2, LAMAX, NAST, RA, BSTO
      if(NAST.ne.ntarg)  Stop 'NAST <> ntarg'

      read(in) (E,N=1,NAST)
      read(in) (L,N=1,NAST)
      read(in) (I,N=1,NAST)
      read(in) ((C,K=1,3),L=1,LRANG2)

! ... Cycle over partial waves

    1 Continue
      read(in,end=2) IL2, IS2, IP2, NCHAN, NHM, MORE
      IP2 = (-1)**IP2

! ... find partial wave:

      klsp = 0
      Do i=1,nlsp
       if(IL2.ne.Lpar(i)) Cycle
       if(IS2.ne.ISpar(i)) Cycle
       if(IP2.ne.ipar(i)) Cycle
       klsp=i; Exit
      End do
      if(klsp.eq.0) Stop 'unknown partial wave in H.DAT'
      if(nch(klsp).ne.nchan) Stop 'nchan <> nch'

      read(in) (NCONAT, N=1,NAST)
      read(in) (L2P, I=1,NCHAN)

      if(allocated(CF)) Deallocate(CF); Allocate (CF(nchan,nchan,lamax))

      read(in) (((CF(I,J,K), I=1,nchan), J=1,nchan), K=1,LAMAX)

      read (in) (W,i=1,NHM)      ! skip eigenvalues
      read (in) ((W,i=1,nchan),j=1,NHM)   ! skip surface amplitudes

! ... calculations of f-values:

      if(klsp1.gt.0.and.klsp.lt.klsp1) go to 1
      if(klsp2.gt.0.and.klsp.gt.klsp2) go to 1

      if(LPAR1.ge.0.and.IL2.lt.LPAR1) go to 1
      if(LPAR2.ge.0.and.IL2.gt.LPAR2) go to 1

      nu = 3
      write(AF,'(a,i3.3)') 'f_values_',klsp
      Open(nu,file=AF)

      eps = 1.d-16
      Call f_values(nut,klsp,np,ni,nu,eps,lamax,nchan,CF,OS)

      if(MORE.gt.0) go to 1

    2 Continue

      AF = 's_values'
      Open(nu,file=AF)
      n = count(OS(:,:).ne.0.d0)
      write(nu,'(i12,a)') n,'    - number of dipole transitions'
      Do i=1,ntarg-1; Do j=i+1,ntarg
       if(OS(i,j).eq.0.d0) Cycle
       write(nu,'(2i6,E15.5)') i,j, OS(i,j)
      End do; End do


      End ! program


!======================================================================
      Subroutine f_values(nut,klsp,np,ni,pri,eps_acf,km,nchan,ACF,OS)
!======================================================================
!     define the f-values between target states based on the
!     given asimptotic coefficients ACF for k=1
!----------------------------------------------------------------------
      Use zconst, only: c_au, time_au
      Use target
      Use channel

      Implicit none
      Real(8) :: ACF(nchan,nchan,km)

      Real(8) :: AK(ntarg,ntarg), RDME(ntarg,ntarg), OS(ntarg,ntarg)
      Integer :: IP(ntarg,ntarg)
      Real(8), external :: Z_6jj, Reduce_factor
      Real(8) :: S,SS, g1,g2, de, a,f, eps_acf
      Integer :: i,j, i1,i2,nt, pri, km, nchan, nut, klsp, k,ij, np,ni
      Character(100) :: line

      Call R_channel(nut,klsp)

      AK=0.d0; IP=0; RDME=0.d0
      Do i=1,nch-1; Do j=i+1,nch
       S=ACF(i,j,1)/2.d0;
       i1=iptar(i); i2=iptar(j)
       if(abs(S).lt.eps_acf) Cycle

       if(i1.gt.ni.and.i2.gt.np) Cycle

       SS = Reduce_factor(i,j,1)

       if(abs(SS).lt.eps_acf) Cycle
       S = S / SS
       de=Etarg(i2)-Etarg(i1)
       RDME(i1,i2) = RDME(i1,i2) + S

       if(istarg(i1).ne.0) then
        S = S*S*istarg(i1)
        g1 = (2*ltarg(i1)+1) * istarg(i1)
        g2 = (2*ltarg(i2)+1) * istarg(i2)
       else
        S = S*S
        g1 = jtarg(i1)
        g2 = jtarg(i2)
       end if

       f = 2.d0/3.d0*de*S /g1
       a = 4.d0/3.d0*de**3*S/c_au**3/time_au /g2

       RDME(i2,i1) = RDME(i2,i1) + S
       AK(i1,i2) = AK(i1,i2) + f
       AK(i2,i1) = AK(i2,i1) + a
       IP(i1,i2) = IP(i1,i2) + 1
      End do; End do

      nt = 0
      Do i1=1,ntarg; Do i2=i1,ntarg
       if(IP(i1,i2).eq.0) Cycle
       RDME(i1,i2) = RDME(i1,i2) / IP(i1,i2)
       RDME(i2,i1) = RDME(i2,i1) / IP(i1,i2)
       AK(i1,i2) = AK(i1,i2) / IP(i1,i2)
       AK(i2,i1) = AK(i2,i1) / IP(i1,i2)
       nt = nt + 1
      End do; End do

! ... total decay probabilities:

      Do i=2,ntarg;  AK(i,i) = SUM(AK(i,1:i-1)); End do

! ... print results:

      if(nt.gt.0) &
      write(pri,'(/a,a/)') 'Target radiative data:  '
      write(pri,'(/a,a/)') 's-value, f-value, A-value, branching ratio,',&
                           ' dipole reduced matrix elements'

      Do i=1,ntarg-1; Do j=i+1,ntarg
       if(AK(i,j).eq.0.d0) Cycle

       write(pri,'(2i5,E15.5, 1PE12.3,E12.3,0Pf10.5,E15.5,5x,a,a)') &
        i,j, RDME(j,i), AK(i,j),AK(j,i),AK(j,i)/AK(j,j), RDME(i,j),BFT(i),BFT(j)

       if(OS(i,j).eq.0.d0) then
        OS(i,j) = RDME(j,i)
       else
        S = abs(OS(i,j) - RDME(j,i)) / abs(OS(i,j) + RDME(j,i))
        if(S.gt.0.01) write(*,'(i5,2i3,2E15.5,f10.2)') klsp,i,j,OS(i,j),RDME(j,i),S*100
       end if

      End do; End do

! ... polarizability of the ground state:

      a = 0.d0
      Do i = 2,ntarg
       de = Etarg(i)-Etarg(1)
       a = a + AK(1,i)/(de*de)
      End do

      if(a.ne.0.d0) &
      write(pri,'(/a,f10.3/)') 'Polarizability of the ground state =',a

      End  Subroutine f_values

!======================================================================
      Real(8) Function Reduce_factor(ich,jch,k)
!======================================================================
!     define factor connecing reduced dipole matrix element with
!     asymptotic coefficient for multipole index k
!----------------------------------------------------------------------
      Use target
      Use channel

      Implicit none
      Integer, intent(in) :: ich,jch,k
      Integer :: ll1,ll2,jj1,jj2,it,jt,L1,L2,S1,S2,J1,J2,LT,JJ, ip,jp,kz
      Real(8) :: S,SS, zero = 0.d0
      Real(8), external :: ZCLKL, Z_6jj, Z_6j

      Reduce_factor = zero

      ll1 = lch(ich);  ll2 = lch(jch)
      S = ZCLKL(ll1,k,ll2)
      if(S.eq.zero) Return

      jj1 = jkch(ich); jj2 = jkch(jch)
      it = iptar(ich); jt = iptar(jch)
      L1 = ltarg(it);  L2 = ltarg(jt)
      S1 = istarg(it); S2 = istarg(jt)
      J1 = jtarg(it);  J2 = jtarg(jt)
      ip = iptarg(it); jp = iptarg(jt)
      kz = 0; if(ip.ne.jp) kz=1;  kz = kz - k
      if(MOD(kz,2).ne.0) then; Reduce_factor = zero; Return; end if

      LT = lpar; JJ =jpar
      if(coupling.eq.'LS') then

       if(S1.ne.S2) then; Reduce_factor = zero; Return; end if
       S = S * Z_6jj(L1,ll1,LT,ll2,L2,k)
       kz = L2+ll1+LT
       S = S * (-1) ** kz

      elseif(coupling.eq.'JK') then

       if(jj1.ne.jj2) Return
       S = S * Z_6j(ll2+ll2+1,ll1+ll1+1,k+k+1,J1,J2,jj1)
       kz = J2+JJ+JJ-ll1-ll1-jj1+1; kz=kz/2
       S = S * (-1) ** kz

      elseif(coupling.eq.'JJ') then

       SS = jj1*jj2;   S = S * sqrt(SS)
       S = S * Z_6j(ll1+ll1+1,2,jj1,jj2,k+k+1,ll2+ll2+1)
       S = S * Z_6j(J1,jj1,JJ,jj2,J2,k+k+1)
       kz = J2+jj1+JJ+ll1+ll1+1+jj2+k+k; kz = kz/2
       S = S * (-1) ** kz

      end if

      Reduce_factor = S

      End Function Reduce_factor



