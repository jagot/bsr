!======================================================================
      PROGRAM test_radial
!======================================================================
!
!   This program tests the radial routines with the logarithmic scale
!   on the set of hydrogenic orbitals
!
!   INPUT:   int_inp    -  required set of integrals
!
!   OUTPUT:  int_out    -  results for intrgrals
!            test_...   -  others test results (see below)
!
!   SUBROUTINE called:
!                          ALOC_RADIAL
!                          INITR
!                          test_integrals
!                          test_quadr
!                          test_hl
!                          test_rlsht
!                          test_deriv
!----------------------------------------------------------------------

      USE RADIAL

      IMPLICIT REAL (8) (A-H,O-Z)

! ... open files:

      in = 1;      Open(in,file='int_inp',status='OLD')
      iout = 2;    Open(iout,file='int_out')

! ... read basis information:

      read(1,*) Z        !  nuclear charge
      read(1,*) n_max    !  max. n for hydrogenic orbitals

! ... allocate arrays:

      krf = (n_max+1)*n_max/2 + 2

      Call ALOC_RADIAL(krf)

! ... set up grid points:

      CALL INITR

! ... define the set of hydrogenic orbitals:

      Call define_orbitals(n_max)

! ... test assignment for fine constant:

      fine = 1.d0

! ... test of integral routines:

      CALL test_integrals(in,iout)

! ... test of QUADR (<.|r^k|.>) and GRAD:

      Call test_quadr

! ... test HL (<.|L|.>):

      Call test_hl

! ... test RLSHT - relativistic shift:

      Call test_rlsht

! ... test of derivitives:

      Call test_deriv


      END PROGRAM test_radial



!======================================================================
      SUBROUTINE define_orbitals(n_max)
!======================================================================
!
!     define set of hydrogenic orbitals for n <= n_max
!
!----------------------------------------------------------------------

      USE RADIAL, L => lro, N => nro, KS => kro, EL => ero, MX => mro

      IMPLICIT REAL (8) (A-H,O-Z)

      Character(4) ELF4, EL4

      nrf = 0
      S = 0.d0
      Do nn=1,n_max
       Do ll=0,nn-1
        nrf = nrf + 1
        i = nrf

        N(i)=nn; L(i)=ll;  KS(i)=0

        EL4 = ELF4(N(i),L(i),KS(i))
        EL(i) = EL4(2:4)

        PN = HNORM(N(I),L(I),Z-S)
        DO J=1,NR-2
         P(j,i) = PN*HWF(N(I),L(I),Z-S,R(J))/R2(J)
         if(abs(P(j,i)).lt.1.d-16) P(j,i)=0.d0
        End do
        P(NR-1:NR,i) = D0

        Do j = NR-2,1,-1
         MX(i) = j
         if(P(j,i).ne.0.d0) Exit
        End do

        AZ(i) = PN*(2.*(Z - 0.5*S)/N(I))**(L(I) + 1)
        mexp(i) = l(i) + 1
        aexp(i) = - D1/mexp(i)
        ZR = Z*R(1)
        bexp(i) = P(1,i)/(AZ(i)*R2(1)*R(1)**L(i))
        bexp(i) = (bexp(i) - D1 + ZR/(L(i)+1) )/ZR**2

       End do
      End do

      END SUBROUTINE define_orbitals



!========================================================================
      SUBROUTINE test_integrals(in,iout)
!========================================================================

! ... tests the given set of radial integrals on the basis of hydrogenic
! ... orbitals. The integrals are given in file 'int_inp', the results are
! ... output in file 'int_out'.

      USE RADIAL

      IMPLICIT REAL *8(A-H,O-Z)

      Character (80) :: AS
      Character (4) :: EL1, EL2, EL3, EL4

      Real(8), External :: NK

    1 read(in,'(a)',end=2) AS

      if(AS(1:1).eq.'*') go to 2          ! end of work

!----------------------------------------------------------------------
      if(AS(1:1).eq.'#') then                                  ! timing

         t1 = RRTC()

         Do k=1,10000
          if(AS(2:3).eq.'RK')   S  = RK  (1,1,1,1,0)
          if(AS(2:3).eq.'NK')   S  = NK  (1,1,1,1,0)
          if(AS(2:3).eq.'VK')   S  = VK  (1,1,1,1,0)
          if(AS(2:3).eq.'TK')   S  = TK  (1,1,1,1,0)
         End do

         t2 = RRTC()

         write(iout,'(/a,a,f10.3,a/)') ' timing ',AS(2:3), &
                                            (t2-t1)/10, '  msec'

         go to 1

      end if
!---------------------------------------------------------------------
!            ... test of accuracy ...

      Read(AS,'(1x,i2,4(1x,a3))') k, el1, el2, el3, el4

      Call EL4_nlk(el1,n1,l1,k1)
      ie1 = IWF(n1,l1,k1)
      Call EL4_nlk(el2,n2,l2,k2)
      ie2 = IWF(n2,l2,k2)
      Call EL4_nlk(el3,n3,l3,k3)
      ie3 = IWF(n3,l3,k3)
      Call EL4_nlk(el4,n4,l4,k4)
      ie4 = IWF(n4,l4,k4)

      read(AS(22:80),*) S1,S2
      SA = S1/S2
      if(AS(1:1).eq.'R') SA = SA * Z
      if(AS(1:1).eq.'N') SA = SA * Z**3
      if(AS(1:1).eq.'V') SA = SA * Z**3
      if(AS(1:1).eq.'T') SA = SA * Z**3 /(2*k+1)

      write(iout,'(a21,d24.16)') AS(1:21),sa

      if(AS(1:1).eq.'R') then

        S  = RK  (ie1,ie2,ie3,ie4,k)
        SS = (S-SA)/SA
        write(iout,'(21x,d24.16,d12.3)') s,ss
        write(*,'(21x,d24.16,d12.3)') s,ss

      elseif(AS(1:1).eq.'N') then

       S  = NK  (ie1,ie2,ie3,ie4,k)
       SS = (S-SA)/SA
       write(iout,'(21x,d24.16,d12.3)') s,ss
       write(*,'(21x,d24.16,d12.3)') s,ss

      elseif(AS(1:1).eq.'V') then

       S  = VK (ie1,ie2,ie3,ie4,k)
       SS = (S-SA)/SA
       write(iout,'(21x,d24.16,d12.3)') s,ss
       write(*,'(21x,d24.16,d12.3)') s,ss

      elseif(AS(1:1).eq.'T') then

       S  = TK (ie1,ie2,ie3,ie4,k)
       SS = (S-SA)/SA
       write(iout,'(21x,d24.16,d12.3)') s,ss
       write(*,'(21x,d24.16,d12.3)') s,ss

      end if

      go to 1

    2 Continue

      END SUBROUTINE test_integrals



!========================================================================
      SUBROUTINE test_quadr
!========================================================================

! ... tests the accuracy of one-electron radial integrals on the basis
! ... of hydrogenic orbitals. Results are output in file 'quadr_out'.

      USE RADIAL

      IMPLICIT REAL *8(A-H,O-Z)

      iout=3;  Open(iout,file='test_quadr')
!----------------------------------------------------------------------
!                                            accuracy of orthogonality:

      write(iout,'(/a/)')   '  accuracy of orthogonality'
      Do i=1,nrf
       Do j=i,nrf

       if(l(i).ne.l(j)) Cycle

       a = QUADR(i,j,0)
       if(i.eq.j) a = a - 1.d0

       write(iout,'(3x,a,3x,a,D12.3)') EL(i),EL(j),a

       End do
      End do
!-----------------------------------------------------------------------
!                                           accuracy of avarages <r^m>:

      write(iout,'(/a/)')   '  accuracy of <r^m>'
      write(iout,'(/4x,6(7x,a5)/)') &
          ' <r> ','<r^2>','<r^3>','<r-1>','<r-2>','<r-3>'

      Do i=1,nrf
       write(iout,*)
       Call HVA(n(i),l(i),z,ar1,ar2,ar3,am1,am2,am3)
       br1 = QUADR(i,i,1)-ar1
       br2 = QUADR(i,i,2)-ar2
       br3 = QUADR(i,i,3)-ar3
       bm1 = QUADR(i,i,-1)-am1
       bm2 = QUADR(i,i,-2)-am2
       bm3 = 0.0
       if(l(i).gt.0) bm3 = QUADR(i,i,-3)-am3
       if(ar1.ne.0.d0) br1=br1/ar1
       if(ar2.ne.0.d0) br2=br2/ar2
       if(ar3.ne.0.d0) br3=br3/ar3
       if(am1.ne.0.d0) bm1=bm1/am1
       if(am2.ne.0.d0) bm2=bm2/am2
       if(am3.ne.0.d0) bm3=bm3/am3
       write(iout,'(a4,1x,6D12.3)') el(i),ar1,ar2,ar3,am1,am2,am3
       write(iout,'(   5x,6D12.3)')       br1,br2,br3,bm1,bm2,bm3
       end do
!-----------------------------------------------------------------------
!                                  accuracy of E1 transition amplitudes:

       write(iout,'(/a/)')  &
              '  accuracy of <i|r|j>:             QUADR  and  GRAD'


       Do icase = 1,3

        if(icase.eq.1)   i = IWF(1,0,0)         ! 1s - np
        if(icase.eq.2)   i = IWF(2,0,0)         ! 2s - np
        if(icase.eq.3)   i = IWF(2,1,0)         ! 2p - ns,nd

       Do j=1,nrf

        if(ABS(l(i)-l(j)).ne.1) Cycle
        if(n(j).le.n(i)) Cycle

        a = n(j)

        if(icase.eq.1) &

        ar1 = 16*a**3.5d0*(a-1)**(a-2.5d0)/(a+1)**(a+2.5d0)

        if(icase.eq.2) &

        ar1 = 2**8*a**3.5d0*sqrt(2*(a**2-1))*(a-2)**(a-3)/(a+2)**(a+3)

        if(icase.eq.3.and.l(j).eq.0) &

        ar1 = 2**8*a**4.5d0/sqrt(6.d0)*(a-2)**(a-3)/(a+2)**(a+3)

        if(icase.eq.3.and.l(j).eq.2) &

        ar1 = 2**10*a**4.5d0*sqrt((a**2-1)/6)* &
                 (a-2)**(a-3.5d0)/(a+2)**(a+3.5d0)

        ar1 =  ar1/z

        br1 =  (QUADR(i,j,1) - ar1)/ar1

        a   = -(1.d0/n(i)**2 - 1.d0/n(j)**2)/2.d0 * z**2

        bm1 =  GRAD(i,j) / a
        bm1 =  (bm1 - ar1)/ar1

        bz1 =  -ZGRAD(i,j) / a
        bz1 =  (bz1 - ar1)/ar1


        write(iout,'(a,a,a,a,a,4D12.3)') &
                    '<',EL(i),'| r |',EL(j),'> = ',ar1,br1,bm1,bz1

       End do
       End do


!-----------------------------------------------------------------------
!                                  accuracy of E2 transition amplitudes:

       write(iout,'(/a/)')  &
              '  accuracy of <i|r^2|j>:             QUADR  and  GRAD2'


       Do icase = 1,3

        if(icase.eq.1)   i = IWF(1,0,0)         ! 1s - np
        if(icase.eq.2)   i = IWF(2,0,0)         ! 2s - np
        if(icase.eq.3)   i = IWF(2,1,0)         ! 2p - ns,nd

       Do j=1,nrf

        if(ABS(l(i)-l(j)).ne.0) Cycle
        if(n(j).le.n(i)) Cycle

        a   = (1.d0/n(i)**2 - 1.d0/n(j)**2)/2.d0 * z**2

        ar2 =  QUADR(i,j,2)
        gr2 =  2*Grad2(i,j)/a
        bm2 =  (gr2 - ar2)/ar2

        write(iout,'(a,a,a,a,a,4D12.3)') &
                    '<',EL(i),'| r^2 |',EL(j),'> = ',ar2,gr2,bm2

       End do
       End do



     END SUBROUTINE test_quadr



!========================================================================
      SUBROUTINE test_hl
!========================================================================
!
!     tests accuracy of one-electron L integrals on the basis
!     of hydrogenic orbitals. Results are output in file 'hl_test'
!
!------------------------------------------------------------------------

      USE RADIAL

      IMPLICIT REAL (8) (A-H,O-Z)

      iout=3;  Open(iout,file='hl_test')

      write(iout,'(/a/)') '   accuracy and symmetry of HL:'

      Do i=1,nrf
       Do j=i,nrf

       if(l(i).ne.l(j)) Cycle

       a=0.d0;  if(i.eq.j) a = z*z/n(i)**2

       b=HL(i,j);  c=HL(j,i)

       if(a.ne.0.d0) b=(b-a)/a
       if(a.ne.0.d0) c=(c-a)/a
       write(iout,'(a,a,a,a,a,3D12.3)')  &
                    '<',EL(i),'| L |',EL(j),'> =',a,b,c

       End do
      End do

      Close(iout)

      END SUBROUTINE test_hl



!======================================================================
      SUBROUTINE test_rlsht
!======================================================================
!
! ... tests the accuracy of relativistic corrections on the basis
! ... of hydrogenic orbitals. Results are output in file 'hl_out'.
!
!----------------------------------------------------------------------

      USE RADIAL

      IMPLICIT REAL (8) (A-H,O-Z)

      iout=3;  Open(iout,file='test_rlsht')


      write(iout,'(/a/)') '   accuracy rel. shift:'


                          !   - 0.5 * (Z/n)^4 * [ 4n/(l+1/2) -3 ]
      Do i=1,nrf

       a = l(i)+0.5
       a = 4*n(i)/a - 3
       a = -0.5 * (Z/n(i))**4 * a


       b = RLSHFT(i,i)    !       * 2.d0

       if(l(i).eq.0) b = b - 0.5*Z*AZ(i)*AZ(i)

       c = (a - b)/a

       write(iout,'(a,a,a,a,a,3D20.8)')  &
                    '<',EL(i),'| rel |',EL(i),'> =',a,b,c

      End do

      Close(iout)

    END SUBROUTINE test_rlsht





!========================================================================
      SUBROUTINE test_deriv
!========================================================================
!
! ... test of p_deriv function
!
!------------------------------------------------------------------------


      USE RADIAL

      IMPLICIT REAL *8(A-H,O-Z)

      iout = 3;   Open(iout,file='test_deriv')

! ... 1s:    P' = -2*r*exp(-r)

      write(iout,'(/a/)') ' 1s:'
      ii = IWF(1,0,0)

      IP = nrf+1
      Call P_derive(II,IP)

      mm = max(mx(ii),mx(ip))
      Do  m = 1,mm
       x = R(m); a = -D2*x*exp(-x)/R2(m); b = P(m,ip)
       write(iout,'(i5,4E15.5)')  m, x,a,b,(a-b)/a
      End do


! ... 2p:    P' =    r*exp(-r/2)(1-r/2)

      write(iout,'(/a/)') ' 2p:'
      ii = IWF(2,1,0)

      IP = nrf+1
      Call P_derive(II,IP)

      mm = max(mx(ii),mx(ip))
      pn = 1.d0/(2.d0*sqrt(6.d0))
      Do  m = 1,mm
       x = R(m); a = pn*x*exp(-x/2)*(1-x/2)/R2(m); b = P(m,ip)
       write(iout,'(i5,4E15.5)')  m, x,a,b,(a-b)/a
      End do

! ... 2s:    P' = -(1/sqrt(2))*r*exp(-r/2)(1-r/4)

      write(iout,'(/a/)') ' 2s:'
      ii = IWF(2,0,0)

      IP = nrf+1
      Call P_derive(II,IP)

      mm = max(mx(ii),mx(ip))
      pn = 1/sqrt(2.d0)
      Do  m = 1,mm
       x = R(m); a = -pn*x*exp(-x/2)*(1-x/4)/R2(m); b = P(m,ip)
       write(iout,'(i5,4E15.5)')  m, x,a,b,(a-b)/a
      End do

      End SUBROUTINE test_deriv





!======================================================================
      Real(8) FUNCTION QUADR(I,J,KK)
!======================================================================
!
!                                 kk
!     Evaluates the integral of  r   P (r) P (r) with respect to r
!                                     i     j
!----------------------------------------------------------------------

      USE RADIAL, MX => mro

      IMPLICIT NONE
      Integer(4), Intent(in) :: i,j,kk
      Integer(4) :: k,m,mm
      REAL(8) :: A,B,DEN,ZR

! ... region (0,r1) ...

      m   = mexp(i) + mexp(j) + kk
      DEN = m + 1
      A  = aexp(i) + aexp(j)
      B  = bexp(i) + bexp(j) + aexp(i)*aexp(j)
      B  = (DEN+D1)*A**2 - D2*B/(DEN+D2)
      A  = - A/(DEN + D1)

      K = KK + 2
      ZR = Z*R(1)
      A  = P(1,I)*P(1,J)*R(1)**K * (((B*ZR + A)*ZR + D1)/(DEN*H1) + D5)

! ... region (r1,inf) ...

      MM = MIN0(MX(I),MX(J))

      DO m = 3,MM,2
       A = A + P(m,I)*P(m,J)*R(m)**K
      END DO

      B = D0
      DO m = 2,MM,2
       B = B + P(m,I)*P(m,J)*R(m)**K
      END DO

      QUADR = H1*(A + B + B)

      END FUNCTION QUADR
