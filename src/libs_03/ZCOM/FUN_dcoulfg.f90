!=====================================================================
      Subroutine DCOULFG(E,K,R,Z,AMASS,FREG,GREG,FIRR,GIRR, &
                         FPREG,GPREG,FPIRR,GPIRR,TAND,MODE)
!=====================================================================
!     Dirac continuum Coulomb Function
!     (after Toffoli and Decleva, CPC, 152 (2003) 151)
!
!     Call: COULFG
!
!     mode = 0 - only function values;
!          = 1 - also derivaties
!     TAND - phase
!---------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)
      Parameter(CI=137.03598950D0, PIG=3.1415926535897930D0)
      Parameter(MLAT=200)
      Dimension FC(mlat),GC(mlat),FCP(mlat),GCP(mlat)

      alpha = 1.d0/CI
      amc2 = amass*CI**2
      akappa = sqrt(E**2-amc2**2)/ci
      AR = R*akappa
      AK = abs(k)
      W = k/AK
      zalpha = Z*alpha
      gamma = sqrt(ak**2-zalpha**2)
      L1 = INT(gamma+1.d0)
      if(L1.gt.mlat) Stop 'DCOULFG: increase MLAT !'

      EFAC = sqrt((E+amc2)/(E-amc2))
      ETA1 = -(E*zalpha)/(CI*akappa)
      A = SQRT((E*AK)**2 - (gamma*amc2)**2)/(E*AK+gamma*amc2)
      B = zalpha/(AK+gamma)
      cnorm = (SQRT((AK+gamma)*(E*AK+gamma*amc2)/(2*PIG*AKAPPA)))/ &
              (gamma*CI)

      SMIN = gamma-1.d0
      M1 = max(INT(SMIN),0)+1
      XMIN = MOD(SMIN,1.d0)
      kfn = 0
      mode1 = 2
      if(mode.eq.1) mode1=1

      Call COULFG(AR,ETA1,XMIN,gamma,FC,GC,FCP,GCP,mode1,kfn,ifail)

      if(ifail.ne.0) Stop 'COULFG: ifail \= 0'

      iup = M1+1
      ilow = M1
      if(k.lt.0) then
       iup=M1
       ilow=M1+1
       EFAC = 1.d0/EFAC
      end if

      FREG = CNORM *(FC(iup)+A*B*FC(ilow))
      FIRR = CNORM *(GC(iup)+A*B*GC(ilow))
      GREG = W*CNORM *(B*FC(iup)+A*FC(ilow))
      GIRR = W*CNORM *(B*GC(iup)+A*GC(ilow))

      TAND = -zalpha * efac /(AK+gamma)

      if(mode.eq.1) then
       FPREG = CNORM *(FCP(iup)+A*B*FCP(ilow))*akappa
       FPIRR = CNORM *(GCP(iup)+A*B*GCP(ilow))*akappa
       GPREG = W*CNORM *(B*FCP(iup)+A*FCP(ilow))*akappa
       GPIRR = W*CNORM *(B*GCP(iup)+A*GCP(ilow))*akappa
      end if

      End Subroutine DCOULFG


!=====================================================================
      Subroutine ZCOUL(E,L,Z,R,F,G,FP,GP)
!=====================================================================
!
!   Gives the magnitude and derivative for continuuum Coulomb Function
!   with k^2 = E, orbital momentum L, and charge Z in point R
!
!   Call: COULFG
!---------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Integer, intent(in) :: L
      Real(8), intent(in) :: E,Z,R
      Real(8), intent(out) :: F,G,FP,GP
      Integer, parameter :: max_lc = 99
      Integer :: lmin, lmax, ifail
      Real(8) :: RHO,ETA,K,AN,FC,FCP,GC,GCP, XL,XM

      COMMON /COUL/ FC(100),GC(100),FCP(100),GCP(100)

      if(l.gt.max_lc) Stop ' ZCFG: l > max_lc=99'

      K = SQRT(E)
      ETA = -Z/K
      RHO = K*R
      XL = L
      XM = L

      CALL COULFG(RHO,ETA,XL,XM,FC,GC,FCP,GCP,1,0,IFAIL)

      AN=DSQRT(1.0/K)
      F = FC(1)*AN
      G = GC(1)*AN
      FP = FCP(1)/AN
      GP = GCP(1)/AN

      End Subroutine ZCOUL
