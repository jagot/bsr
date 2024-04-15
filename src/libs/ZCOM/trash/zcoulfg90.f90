!======================================================================
      SUBROUTINE ZCOULFG90 (X,ETA,XL, FC,GC,FCP,GCP, KFN, ifAIL)
!======================================================================
!
!     COULOMB & BESSEL FUNCTION PROGRAM USING STEED'S METHOD
!
!   COUL90 RETURNS ARRAYS FC = F, GC = G, FCP = F'(x), GCP = G'(x)
!   FOR REAL X .GT. 0., REAL ETA (INCLUDING 0.), AND REAL XL .GT.-1.
!   IT HENCE GIVES POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER
!   EQUATION, TO THE KLEIN-GORDoN EQUATION AND TO SUITABLE FORMS OF
!   THE DIRAC EQUATION.    BY SETTING ETA = 0.0 AND RENORMALISING
!   SPHERICAL & CYLINDRICAL BESSEL FUNCTIONS ARE COMPUTED INSTEAD.
!  --------------------------------------------------------------------
!   CALLING VARIABLES:
!
!   X       - REAL ARGUMENT FOR COULOMB FUNCTIONS > 0.0 ( x = kr)
!             [ X > SQRT(ACCUR): ACCUR IS TARGET ACCURACY = 1.0D-14 ]
!   ETA     - REAL SOMMERFELD PARAMETER z/k, UNRESTRICTED, > = < 0.0
!   XL      - REAL LAMBDA-VALUE (L-VALUE OR ORDER),
!   FC ,GC  - REAL VALUES F,G OF REGULAR, IRREGULAR COULOMB FUNCTIONS
!   FCP,GCP - REAL VECTORS FOR THE X-DERIVATIVES OF  F,G
!             THESE VECTORS TO BE OF LENGTH AT LEAST MINL + LRANGE
!             STARTING ELEMENT MINL = MAX0( IDINT(XLMIN+ACCUR),0 )
!   KFN     - Integer CHOICE OF FUNCTIONS TO BE COMPUTED:
!           = 0         REAL COULOMB FUNCTIONS AND DERIVATIVES F & G
!           = 1    SPHERICAL BESSEL      "      "     "        j & y
!           = 2  CYLINDRICAL BESSEL      "      "     "        J & Y
!
!   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
!   IN OSCILLATING REGION X .GE. [ETA + SQRT{ETA**2 + XLM*(XLM+1)}]
!   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100 * ACC WHERE ACC IS
!   THE SMALLEST NUMBER WITH 1.+ACC.NE.1.
!   if X IS SMALLER THAN [ ] ABOVE THE ACCURACY BECOMES STEADILY WORSE
!  --------------------------------------------------------------------
!   ERROR RETURNS                THE USER SHOULD TEST ifAIL ON EXIT
!
!   ifAIL ON INPUT IS SET TO 0                        LIMIT = 20000
!   ifAIL IN OUTPUT =  0 : CALCULATIONS SATISFACTORY
!                   =  1 : F'/F DID NOT CONVERGE AFTER LIMIT ITERATIONS
!                   =  2 : CF2 DID NOT CONVERGE AFTER LIMIT ITERATIONS
!                   = -1 : X < SQRT(ACCUR)
!                   = -2 : INCONSISTENCY IN ORDER VALUES (L-VALUES)
!  --------------------------------------------------------------------
!  MACHINE-DEPENDENT PARAMETERS:    ACCUR - SEE ABOVE
!           SMALL - OFFSET FOR RECURSION = APPROX SQRT(MIN REAL NO.)
!  --------------------------------------------------------------------
!  PROGRAMMING HISTORY AND BACKGROUND: CPC IS COMPUTER PHYSICS COMMUN.
!  ORIGINAL PROGRAM  RCWFN       IN    CP!  8 (1974) 377-395
!                 +  RCWFF       IN    CPC 11 (1976) 141-142
!  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314
!  REVISED STANDARD  COULFG      IN    CPC 27 (1982) 147-166
!  BACKGROUND MATERIAL IN J. COMP. PHYSICS 46 (1982) 171-188
!  CURRENT PROGRAM   COUL90  (FORTRAN77) SUPERCEDES THESE EARLIER ONES
!  (WHICH WERE WRITTEN IN FORTRAN 4) AND ALSO BY INCORPORATING THE NEW
!  LENTZ-THOMPSON ALGORITHM FOR EVALUATING THE FIRST CONTINUED FRACTION
!  ..SEE ALSO APPENDIX TO J. COMP. PHYSICS 64 (1986) 490-509     1.4.94
!  .. K,Barschat(ed.) Computational Atomic Physics (Springer,1996)
!  --------------------------------------------------------------------
!  AUTHOR: A. R. BARNETT           MANCHESTER  MARCH   1981/95
!                                  AUCKLAND    MARCH   1991
!  ADOPT by Oleg Zatsarinny
!  --------------------------------------------------------------------
!  CALLS:  DSQRT,DABS,MAX0,IDINT,DSIGN,DFLOAT,DMIN1
!  --------------------------------------------------------------------

      Implicit none
      Integer, intent(in)  :: KFN
      Integer, intent(out) :: ifAIL
      Real(8), intent(in)  :: X, ETA, XL
      Real(8), intent(out) :: FC, GC, FCP, GCP
      Real(8) :: ACCUR, ACCH, SMALL,XLM, F,DCF, FCL,GCL,FJWKB,GJWKB, &
                 XI,ETAL,ETAK,RN,TK,DEN,PK,QK, P,Q,                  &
                 WI, A,B,C,D, AR,AI,BR,BI,DR,DI,DP,DQ,               &
                 ALPHA,BETA,GAMMA,OMEGA
      Integer :: IEXP, I
      Logical :: XLTURN
      Real(8), parameter ::  ZERO = 0.d0, ONE = 1.d0, TWO = 2.d0, &
                             TEN2 = 100.d0, HALF = 0.5d0
      Integer, parameter :: LIMIT = 20000

      ACCUR  = TEN2*EPSILON(one)     !  ~ 10^-14
      ACCH   = DSQRT(ACCUR)          !  ~ 10^-7
      SMALL  = DSQRT(TINY(one))      !  ~ 10^-150

      ifAIL = 0
      IEXP  = 1
      GJWKB = ZERO

! ... Check input data: ETA

      if(KFN .NE. 0 .and. ETA .NE. ZERO) then
        WRITE (6,'(/a/)') ' COUL90: ETA <> 0  for kfn > 0'
        ifail = -3
        RETURN
      end if

! ... TEST RANGE OF X, EXIT if.LE.DSQRT(ACCUR) OR if NEGATIVE

      if( X .LE. ACCH ) then
        ifAIL = -1
        WRITE(6,'(/a/)') ' COUL90: X is too small or negative '
        WRITE(6,'(a,D12.3)') ' X = ',X
        WRITE(6,'(a,D12.3)') ' SQUARE ROOT(ACCURACY) =  ',ACCH
        RETURN
      end if

! ... Check input data: XL

      XLM = XL
      if( KFN.EQ.2 ) XLM = XL - HALF
      if( XLM.LE.-ONE) then
        ifAIL = -2
        WRITE (6,'(/a/)') ' COUL90: PROBLEM WITH INPUT L VALUE'
        WRITE (6,'(a,I10,2D15.6)')' XLM = ', XLM
        RETURN
      end if

! ... EVALUATE   F'(L,ETA,X) / F(L,ETA,X)  --> F

      ETAL   = XLM * XLM + XLM
      XLTURN = X * (X -  TWO * ETA) .LT. ETAL
      ETAL   = ETAL  +  ETA * ETA

      XI   = ONE / X
      DEN  = ONE                            ! unnormalised F(XL,ETA,X)
      PK   = XLM + ONE
      F    = ETA / PK  +  PK * XI
      if(DABS(F).LT.SMALL)  F = SMALL
      RN   = ONE
      D    = ZERO
      C    = F

! ... BEGIN LENTZ-THOMPSON PROCEDURE FOR F'/F:

      Do I = 1,LIMIT
        QK = PK + ONE
        if( ETA .NE. ZERO  ) then
          ETAK = ETA / PK
          RN   = ONE + ETAK * ETAK
          TK   = (PK + QK) * (XI + ETAK / QK)
        else
          TK   = (PK + QK) * XI
        end if
        D   =  TK - RN * D          ! direct  ratio of B convergents
        C   =  TK - RN / C          ! inverse ratio of A convergents
        if( DABS(C).LT.SMALL ) C = SMALL
        if( DABS(D).LT.SMALL ) D = SMALL
        D   = ONE / D
        DCF = D * C
        F   = F * DCF
        if( D.LT.ZERO )    DEN = - DEN
        PK  = QK
        if( DABS(DCF-ONE).LT.ACCUR )  EXIT
      End Do

      if( DABS(DCF-ONE).ge.ACCUR ) then
       ifAIL =  1
       WRITE(6,'(/a,I10,a/)') ' COUL90: F''/F HAS FAILED TO CONVERGE AFTER',&
                                LIMIT,'  ITERATIONS '
       WRITE(6,'(a,4D12.3/)') ' F,DCF,PK,ACCUR =  ', &
                                F,DCF,PK,ACCUR
       RETURN
      end if

! ...  EVALUATE P + I.Q  USING STEED'S ALGORITHM (NO ZEROS)

      if( XLTURN ) CALL JWKB90( X,ETA,DMAX1(XLM,ZERO),FJWKB,GJWKB,IEXP )

      if( IEXP.GT.1 .OR. GJWKB.GT.(ONE / (ACCH*TEN2)) ) then
          OMEGA = FJWKB
          GAMMA = GJWKB * OMEGA
          P     = F
          Q     = ONE
      else
          XLTURN = .FALSE.
          PK =  ZERO
          WI =  ETA + ETA
          P  =  ZERO
          Q  =  ONE - ETA * XI
          AR = -ETAL
          AI =  ETA
          BR =  TWO * (X - ETA)
          BI =  TWO
          DR =  BR / (BR * BR + BI * BI)
          DI = -BI / (BR * BR + BI * BI)
          DP = -XI * (AR * DI + AI * DR)
          DQ =  XI * (AR * DR - AI * DI)

          Do I = 1, LIMIT
             P  = P  + DP
             Q  = Q  + DQ
             PK = PK + TWO
             AR = AR + PK
             AI = AI + WI
             BI = BI + TWO
             D  = AR * DR - AI * DI + BR
             DI = AI * DR + AR * DI + BI
             C  = ONE / (D * D + DI * DI)
             DR =  C * D
             DI = -C * DI
             A  = BR * DR - BI * DI - ONE
             B  = BI * DR + BR * DI
             C  = DP * A  - DQ * B
             DQ = DP * B  + DQ * A
             DP = C
             if( DABS(DP)+DABS(DQ).LT.(DABS(P)+DABS(Q)) * ACCUR ) EXIT
          End Do

       if( DABS(DP)+DABS(DQ).GE.(DABS(P)+DABS(Q)) * ACCUR ) then
        ifAIL =  2
        WRITE(6,'(/a,i7,a/)') ' COUL90: CF2 HAS FAILED TO CONVERGE AFTER ',&
                                LIMIT,' ITERATIONS'
        WRITE(6,'(a,1P4D17.7,D12.3/)')' P,Q,DP,DQ,ACCUR =  ', &
                                        P,Q,DP,DQ,ACCUR
        RETURN
       end if


! ... SOLVE FOR FCL = F AT LAMBDA = XLM AND NORMALISING FACTOR OMEGA

      GAMMA   = (F - P) / Q
      if( DABS(GAMMA) .LE. ONE )  then
         OMEGA  = DSQRT( ONE  +  GAMMA * GAMMA )
      else
         OMEGA  = DSQRT( ONE  +  ONE/(GAMMA*GAMMA)) * DABS(GAMMA)
      end if

      OMEGA  = ONE / ( OMEGA * DSQRT(Q) )

     end if

! ... RENORMALISE if SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS

      if( KFN.EQ.1 ) then               ! spherical Bessel functions
        ALPHA = XI
        BETA  = XI
      else if( KFN.EQ.2 ) then          ! cylindrical Bessel functions
        ALPHA = HALF * XI
        BETA  = DSQRT( XI / ASIN(one))! sqrt(2/pi)
      else
        ALPHA = ZERO                    ! kfn = 0,   Coulomb functions
        BETA  = ONE
      end if

      FCL = DSIGN(OMEGA,DEN) * BETA

      if( XLTURN )   then
        GCL =  GJWKB * BETA
      else
        GCL =  FCL * GAMMA
      end if

      if( KFN.NE.0 )    GCL = - GCL     ! Bessel sign differs

      FC  = FCL
      GC  = GCL
      GCP = GCL * (P - Q/GAMMA - ALPHA)
      FCP = FCL * (F - ALPHA)

      END SUBROUTINE ZCOULFG90

!======================================================================
      SUBROUTINE  JWKB90 (X,ETA,XL, FJWKB,GJWKB, IEXP)
!======================================================================
!   COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS  FOR XL .GE. 0.
!   AS MODifIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
!   AUTHOR:    A.R.BARNETT   FEB 1981    LAST UPDATE MARCH 1991
!
!   CALL: DEXP, DLOG, IDINT, MAXEXPONENT, DATAN2, DMAX1, DSQRT, DFLOAT
!  --------------------------------------------------------------------
      Implicit none
      Real(8), intent(in)  :: X,ETA,XL
      Real(8), intent(out) :: FJWKB,GJWKB
      Integer, intent(out) :: IEXP
      Real(8) :: GHH,XLL,HLL,HL,SL,RL,GH,PHI,PHI10
      Real(8), parameter :: ZERO = 0.d0, HALF = 0.5d0, RL35 = 35.d0, &
                            ONE = 1.d0, SIX = 6.d0, TEN = 10.d0

      GHH   =  X * (ETA + ETA - X)
      XLL   = DMAX1( XL * XL + XL, ZERO )

      if( GHH + XLL .LE. ZERO )  RETURN

      HLL  = XLL + SIX / RL35
      HL   = DSQRT(HLL)
      SL   = ETA / HL + HL / X
      RL   = ONE + ETA * ETA / HLL
      GH   = DSQRT(GHH + HLL) / X

      PHI  = X*GH - HALF*( HL*DLOG((GH + SL)**2 / RL) - DLOG(GH) )

      if ( ETA.NE.ZERO ) PHI = PHI - ETA * DATAN2(X*GH, X-ETA)

      PHI10 = - PHI * DLOG(EXP(ONE))
      IEXP  =  IDINT(PHI10)

      if ( IEXP.GT.MAXEXPONENT(one)) then
           GJWKB = TEN**(PHI10 - DFLOAT(IEXP))
      else
           GJWKB = DEXP(-PHI)
           IEXP  = 0
      end if

      FJWKB = HALF / (GH * GJWKB)

      END SUBROUTINE JWKB90



