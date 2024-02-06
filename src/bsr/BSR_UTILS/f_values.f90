!=======================================================================
!     Program f-values
!=======================================================================
! ... calculation of f-values for transition between target states
! ... based on the asymptotic coefficients in the H.DAT file
!-----------------------------------------------------------------------
      Use target
      Use channel

      Implicit real(8) (A-H,O-Z)
      Character(80) :: AF
      Real(8), allocatable :: CF(:,:,:)

!----------------------------------------------------------------------
      iarg = COMMAND_ARGUMENT_COUNT()
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)

      if(AF.eq.'?') then
        write(*,'(a)')
        write(*,'(a)') 'calculation of f-values for transition between target states'
        write(*,'(a)') 'based on the asymptotic coefficients in the H.DAT file'
        write(*,'(a)')
        write(*,'(a)') 'h.nnn or H.DAT + target -->  fvalues.nnn '
        write(*,'(a)')
        write(*,'(a)') 'Call as:    f_values [h=.. klsp=..]'
        write(*,'(a)')
        write(*,'(a)') 'h [H.DAT]  -  name of specific H.DAT-file'
        write(*,'(a)') 'klsp       -  h = h.nnn,  nnn=klsp'
        write(*,'(a)')
        Stop
      end if
!----------------------------------------------------------------------

      klsp=1; Call Read_iarg('klsp',klsp)

! ... target information:

      nut = 1
      Open(nut,file='target',status='OLD')
      Call R_target(nut)
      Call R_channel(nut,klsp)
      close(nut)

! ... read  assymtotic coefficients from H.DAT file

      in = 2
      write(AF,'(a,i3.3)') 'h.',klsp
      if(Icheck_file(AF).eq.0) AF = 'H.DAT'
      Call Read_aarg('h',AF)
      Call Check_file(AF)
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
      Use target
      Use channel

      Implicit none
      Real(8) :: ACF(nch,nch,km)

      Real(8) :: AK(ntarg,ntarg), RDME(ntarg,ntarg)
      Integer :: IP(ntarg,ntarg)
      Real(8), external :: CLEBCH, Z_6jj, Reduce_factor
      Real(8) :: S,SS, g1,g2, de, a,f, eps_acf
      Integer :: i,j, i1,i2,nt, pri, km

      AK=0.d0; IP=0; RDME=0.d0
      Do i=1,nch-1; Do j=i+1,nch; if(i.eq.j) Cycle
       S=ACF(i,j,1)/2.d0;
       i1=iptar(i); i2=iptar(j)
       if(abs(S).lt.eps_acf) Cycle
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
       RDME(i2,i1) = RDME(i2,i1) / IP(i1,i2)
       RDME(i1,i2) = RDME(i1,i2) / IP(i1,i2)
       AK(i1,i2) = AK(i1,i2) / IP(i1,i2)
       AK(i2,i1) = AK(i2,i1) / IP(i1,i2)
       nt = nt + 1
      End do; End do

! ... total decay probabilities:

      Do i=2,ntarg;  AK(i,i) = SUM(AK(i,1:i-1)); End do

! ... print results:

      if(nt.gt.0) &
      write(pri,'(/a,i6)') 'Target radiative data:  nt = ',nt
      write(pri,'(/a/)') 'indexes, states, weights, s-value, f-value, A-value:'

      Do i=1,ntarg-1; Do j=i+1,ntarg
       if(AK(i,j).eq.0.d0) Cycle

       if(istarg(i).ne.0) then
        i1 = (2*ltarg(i)+1) * istarg(i)
        i2 = (2*ltarg(j)+1) * istarg(j)
       else
        i1 = jtarg(i)
        i2 = jtarg(j)
       end if

       write(pri,'(2i5,E15.5, 1PE12.3,E12.3,0Pf10.5,E15.5,5x,2a20)') &
        i,j, RDME(j,i), AK(i,j),AK(j,i),AK(j,i)/AK(j,j), RDME(i,j), &
        trim(BFT(i)),trim(BFT(j))

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
      Integer :: ll1,ll2,jj1,jj2,it,jt,L1,L2,S1,S2,J1,J2,LT,JJ, kz
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

      LT = lpar; JJ =jpar

      if(coupling.eq.'LS') then

       if(S1.ne.S2) Return
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
