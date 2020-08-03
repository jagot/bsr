!--------------------------------------------------------------------
      Subroutine SBESJ (X,LMAX,IFAIL,ACCUR,XJ)
!--------------------------------------------------------------------
!     REGULAR SPHERISAL BESSEL FUNSTIONS FROM L=0 TO L=LMAX
!     SPECIAL CASE OF COULOMB FUNSTIONS (ETA=0) SEE BARNETT CPC 1980
!--------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Dimension XJ(LMAX+1)
      Data ZERO, ONE /0.0D0,1.0D0/

      IFAIL= 0
      ACC  = ACCUR
      if(DABS(X).LT. ACCUR) then
       IFAIL  = 1
       write(*,'(a,D15.5,a)') 'SBESJ FAILED: X = ',X
       write(*,'(a)') 'TRY ASYMPTOTIC SOLUTION'
       Stop
      end if

      XI   = ONE/X
      W    = XI + XI
      F    = ONE
      FP   = (LMAX + 1)*XI
      B    = FP + FP + XI
      D    = ONE/B
      DEL  = -D
      FP   = FP + DEL
      END  = B + 50000.0*W      ! 2  ???
      Do
       B    = B + W
       D    = ONE/(B - D)
       DEL  = DEL*(B*D - ONE)
       FP   = FP + DEL
       if(D.LT.ZERO) F = -F
       if(B.GT.END) then
        IFAIL  = 1
        write(*,'(a,D15.5,a)') 'SBESJ FAILED: X = ',X
        write(*,'(a)') 'TRY ASYMPTOTIC SOLUTION'
        Stop
       end if 
       if(DABS(DEL).LT.DABS(FP)*ACC) Exit
      End do
      FP   = FP*F
      if(LMAX.EQ.0) GO TO 3
      XJ(LMAX+1) = F
      XP2        = FP

! ... DOWNWARD RECURSION TO L=0 (N.B. COULOMB FUNCTIONS)

      PL   = LMAX*XI
      L    = LMAX
      Do LP = 1,LMAX
       XJ(L)= PL*XJ(L+1) + XP2
       FP   = PL*XJ(L)   - XJ(L+1)
       XP2  = FP
       PL   = PL - XI
       L    =  L - 1
      End do
      F = XJ(1)

! ... SOLVE FOR THE L=0 COULOMB FUNCTIONS

    3 W = XI/DSQRT(FP*FP + F*F)
      XJ(1)= W*F
      Do L = 1,LMAX; XJ(L+1)= W*XJ(L+1); End do

      End Subroutine SBESJ


!====================================================================
      Subroutine SPBES (X,LMAX,F)
!====================================================================
!     (present by A.Grum-Grzhimajlo, Moscow university)
!--------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Dimension F(LMAX)
      Real(8), parameter :: acc = 1.d-30

      LMAX1=(X+0.4d0)/0.8d0
      if(X.GT.40.d0.AND.LMAX1.LT.LMAX) GO TO 2

      LMAX1=MIN0(LMAX,LMAX1); CL=1.d0; L1=1

      if(X.GT.1.d0) then
       FL=SIN(X)/X; F1=COS(X)/X
       Do L=1,LMAX1
        F(L)=FL; F2=F1; F1=FL
        FL=F1*(2.*REAL(L)-1.)/X-F2
       End do
       IF(LMAX1.GE.LMAX) RETURN
       L1=LMAX1+1
       Do L=1,LMAX1; CL=CL*X/(2.*L+1.); End do
      end if

      XX=X*X
      Do L=L1,LMAX
       AL=2*L; AS=1.d0; FL=1.d0; S=0.d0
       Do 
        S=S+1.d0; AS=-AS*XX/(2*S*(AL+2*S-1)); FL=FL+AS
        if(abs(AS).lt.acc) Exit
       End do
       F(L)=FL*CL; CL=CL*X/(AL+1.)
      End do
      Return

2     LMAX1=(X+.37)/.67;  LMAX1=MIN0(LMAX,LMAX1)
      DX=X
      DFL=DSIN(DX)/DX
      DF1=DCOS(DX)/DX
      Do L=1,LMAX1
       F(L)=DFL; DF2=DF1; DF1=DFL; DFL=(2*L-1)/DX*DF1-DF2
      End do
      if(LMAX1.GE.LMAX) Return

      L1=LMAX1+1; DXX=DX*DX; DCL=1.D0
      Do L=1,LMAX1; DCL=DCL*DX/(2.D0*REAL(L)+1.D0); End do

      Do L=L1,LMAX
       AL=2*L; DAS=1.d0; DFL=1.d0; S=0.d0
       Do
        S=S+1; CS=2*S*(AL+2*S-1); DAS=-DAS*DXX/DBLE(CS)
        DFL=DFL+DAS
        if(abs(DAS).lt.acc) Exit
       End do
       F(L)=DFL*DCL; DCL=DCL*DX/(AL+1)
      End do

      End Subroutine SPBES 


