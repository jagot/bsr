       REAL*8 FUNCTION F21(A,B,C,D,EPS,IFAIL)
       IMPLICIT REAL*8 (A-H,O-Z)
C
       IPRINT=IFAIL
       IFAIL=0
       T=(A*B*D)/C
       DD=1.0D0/(1.0D0-D)
       SUM=1.0D0+T
       TN1=0.0D0
       I=1
    3  AI=I
       T=T*(A+AI)*(B+AI)*D/((C+AI)*(1.0D0+AI))
       TN2=T*DD
       F21=SUM+TN2
       SUM=SUM+T
       AT=ABS(T+TN2-TN1)
       AS=ABS(F21)*EPS
       IF(AS-AT)1,2,2
    1  TN1=TN2
       I=I+1
       IF(I-300)3,3,4
    4  IF(IPRINT.GT.0)WRITE(6,100)
       IFAIL=3
  100  FORMAT(' FAILED TO CONVERGE IN F21')
    2  RETURN
      END
C
C**********************************************************
C
      REAL*8 FUNCTION FDIP(EK1,L1,EK2,L2,IFAIL,TZED)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C ALAN BURGESS DEPT. OF APPLIED MATHS. AND THEORETICAL PHYSICS,CAMBRIDGE
C CALCULATES THE FUNCTION I(KAPPA1,L1,KAPPA2,L2,1) DEFINED IN PHIL.
C TRANS. ROY. SOC. A226,255,1970, WHERE EK1=KAPPA1**2 AND EK2=KAPPA2**2.
C IT IS SUITABLE FOR USE IN EQUATIONS (8),(9),(10) OR (11) OF
C J. PHYS. B. 7,L364,1974.
C NRB - IFAIL
C      COMMON/NRBZED/TZED,LPRTSW
      DATA EPS/1.D-4/
C
      IF(TZED.EQ.0)THEN
        FDIP=FDIP0(EK1,L1,EK2,L2,EPS,IFAIL)
        RETURN
      ENDIF
C
      IPRINT=IFAIL
      IF(EK1+EK2-1.0D-40) 11,11,12
   11   FDIP=0.0D0
        IFAIL=1
        IF(IPRINT.GT.0)WRITE(6,100)IFAIL
        RETURN
   12 IF(EK1-EK2) 1,1,2
    1 EMIN=EK1
      EMAX=EK2
      GO TO 3
    2 EMIN=EK2
      EMAX=EK1
    3 T=EMIN/EMAX
      IF(T-0.02944D0) 4,4,5
    4 FDIP=FDIP1(EK1,L1,EK2,L2)
      GO TO 9
    5 IF(T-0.16667D0) 7,6,6
    6 FDIP=FDIP2(EK1,L1,EK2,L2)
      GO TO 9
    7 FDIP=FDIP1(EK1,L1,EK2,L2)
      IF(FDIP*FDIP-1.0D-40) 6,6,8
    8   IF(FDIP.LT.0.0.OR.FDIP.GT.1.)THEN
          IFAIL=3
          IF(IPRINT.GT.0)WRITE(6,100)IFAIL
          FDIP=0.0
          RETURN
        ENDIF
        FA=FDIPA(EK1,L1,EK2,L2)
        IFAIL=0
        IF(FA.EQ.0.0)THEN
          FA=FDIP0(EK1,L1,EK2,L2,EPS,IFAIL)
          IFAIL=-IFAIL
          IF(FA.EQ.0.0)RETURN
        ENDIF
        RAT=FDIP/FA
        IF(RAT.GT.10.)THEN
          IFAIL=4
          IF(IPRINT.GT.0)WRITE(6,100)IFAIL
          FDIP=0.0
        ENDIF
        RETURN
    9 IF(FDIP*FDIP-1.0D-40) 10,10,8
   10   IFAIL=2
        IF(IPRINT.GT.0)WRITE(6,100)IFAIL
        RETURN
  100 FORMAT('***FDIP FAILURE: IFAIL=',I2)
      END
C***********************************************************************
      REAL*8 FUNCTION FDIP0(EK1,L1,EK2,L2,EPS,IFAIL)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C ALAN BURGESS,DEPT OF APPLIED MATHS. AND THEORETICAL PHYSICS,CAMBRIDGE
C CALCULATES THE FUNCTION I0(K1,L1,K2,L2,1) DEFINED IN PHIL. TRANS.
C ROY. SOC. A266,255,1970, WHERE EK1=K1*K1, EK2=K2*K2, AND THE RELATIVE
C ACCURACY IS APPROXIMATELY EPS.
C IT IS SUITABLE FOR USE IN EQUATIONS (13) ETC. OF J.PHYS.B. 7,L364,1974
C NRB - IFAIL
C
      IPRINT=IFAIL
      IFAIL=0
      IF(L1-L2)1,2,4
1     L=L1
      GO TO 5
2     IF(IPRINT.GT.0)WRITE(6,100)L1
      IFAIL=1
      FDIP0=0.0D0
3     RETURN
4     L=L2
5     EL=L
      FDIP0=0.5D0/(EL+1.0D0)
      IF(EK1-EK2)6,3,7
6     E=EK1/EK2
      P=L1-L
      GO TO 8
7     E=EK2/EK1
      P=L2-L
8     FDIP0=FDIP0*E**((EL+P+0.5D0)*0.5D0)
C TO OBTAIN THE FUNCTION EK1 OF M.J. SEATON, PROC. PHYS. SOC. A68,457,
C 1955, REMOVE THE 'C' ON THE NEXT LINE.
C     FDIP0=1.0D0
      IF(E -0.5D0)21,20,20
20    P1=P-0.5D0
      T=P1*(EL+1.0D0)*(E -1.0D0)
      I0=L+1
      H0=0.0D0
      DO 9 I=1,I0
      TI=I
      H0=H0+1.0D0/TI
9     CONTINUE
      X=1.0D0-E
      H=1.0D0-(P+P+H0+LOG(0.25D0*X))
      S=1.0D0+T*H
      A=EL+1.0D0
      B=P1
      C=1.0D0
      D=0.0D0
10    A=A+1.0D0
      B=B+1.0D0
      C=C+1.0D0
      D=D+1.0D0
      T=T*A*B*X/(C*D)
      H=H+P1/(D*B)+EL/(C*A)
      T1=T*H
      S=S+T1
      IF(ABS(T1)-EPS*ABS(S))13,11,11
11    IF(C-300.0D0)10,12,12
12    IF(IPRINT.GT.0)WRITE(6,101)
      IFAIL=2
13    FDIP0=FDIP0*S
      RETURN
21    A=EL+1.0D0
      B=P-0.5D0
      C=EL+P+1.5D0
      F=F21(A,B,C,E,EPS,IFAIL)
      L=L+1
      EL=L
      IF(P-0.5D0)23,23,24
23    C1=EL+EL+1.0D0
      GO TO 25
24    C1=1.0D0
25    DO 22 I=1,L
      AI=I
      AII=AI+AI
      C1=C1*AI*AI*4.0D0/(AII*(AII+1.0D0))
22    CONTINUE
      FDIP0=FDIP0*F*C1
      RETURN
100   FORMAT(' FAILED IN FDIP0, L1=L2=',I5)
101   FORMAT(' FAILED TO CONVERGE IN FDIP0')
      END
C***********************************************************************
       REAL*8 FUNCTION FDIP1(EK1,L1,EK2,L2)
C
       IMPLICIT REAL*8 (A-H,O-Z)
C
       IF(L1-L2)1,2,3
    1  L=L1
       A1=EK1
       A2=EK2
       GO TO 4
    2  FDIP1=0.0D0
       RETURN
    3  L=L2
       A1=EK2
       A2=EK1
    4  LP=L+1
       ELP=LP
       B1=SQRT(1.0D0+ELP*ELP*A2)*FMON1(EK1,EK2,L)
       B2=SQRT(1.0D0+ELP*ELP*A1)*FMON1(EK1,EK2,LP)
       IF(B1*B2-1.0D-40)5,5,6
    5  FDIP1=0.0D0
       RETURN
    6  FDIP1=(B1-B2)/ELP
       RETURN
      END
C***********************************************************************
       REAL*8 FUNCTION FDIP2(EK1,L1,EK2,L2)
C
       IMPLICIT REAL*8 (A-H,O-Z)
C
       WMAX=200.0D0
       ETA1=1.0D0/SQRT(EK1)
       ETA2=1.0D0/SQRT(EK2)
       W1=ETA2-ETA1
       PI=3.141592653589793D0
       A=ABS(W1)
       B=PI*A
       IF(B-0.01D0)1,1,2
    1  C=3.0D0/(3.0D0-B*(3.0D0-B*(2.0D0-B)))
       C=SQRT(C)
       GO TO 5
    2  IF(B-14.0D0)4,3,3
    3  C=SQRT(B+B)
       GO TO 5
    4  B=B+B
       C1=1.0D0-EXP(-B)
       C=SQRT(B/C1)
    5  C=0.5D0*C/SQRT(ETA1*ETA2)
       C2=ETA1+ETA2
       C1=4.0D0*ETA1*ETA2/(C2*C2)
       L=L1
       IF(L2-L1)6,6,7
    6  L=L2
       T1=ETA1
       ETA1=ETA2
       ETA2=T1
       W1=-W1
    7  C=C*C1**(L+1)
       U0=L+1
       U1=ETA1
       V0=U0
       V1=-ETA2
       W0=1.0D0
       X0=W1/(C2*C2)
       Y2=-ETA2-ETA2
       Y0=-U0*W1+Y2
       Y1=ETA2*W1
       T1=X0/(1.0D0+W1*W1)
       Z0=U0*T1
       Z1=U1*T1
       T=Z0-Z1*W1
       Z1=Z0*W1+Z1
       Z0=T
       Q0=-1.0D0+Z0*Y0-Z1*Y1
       Q1=Z0*Y1+Z1*Y0
       X=W1*X0
    8  U0=U0+1.0D0
       V0=V0+1.0D0
       W0=W0+1.0D0
       IF(W0-WMAX)21,21,20
   20  FDIP2=0.0D0
       RETURN
   21  CONTINUE
       Y0=Y0+Y2
       T=Z0*U0-Z1*U1
       Z1=Z0*U1+Z1*U0
       Z0=T
       T=Z0*V0-Z1*V1
       Z1=Z0*V1+Z1*V0
       Z0=T
       T=Z0*W0-Z1*W1
       Z1=Z0*W1+Z1*W0
       Z0=T
       X0=X/(W0*(W0*W0+W1*W1))
       Z0=Z0*X0
       Z1=Z1*X0
       T0=Z0*Y0-Z1*Y1
       T1=Z0*Y1+Z1*Y0
       Q0=Q0+T0
       Q1=Q1+T1
       T1=T0*T0+T1*T1
       T0=Q0*Q0+Q1*Q1
       IF(T0-1.0D+24*T1)8,8,9
    9  J1=0
       J2=L+1
       P=ARGAM(J1,W1)+ARGAM(L,ETA1)-ARGAM(J2,ETA2)
       IW0=W0
       IF(A-1.0D-40)11,11,10
   10  P=P+W1*LOG(C2/A)
   11  P0=COS(P)
       P1=SIN(P)
       T=P0*Q0-P1*Q1
       Q1=P0*Q1+P1*Q0
       Q0=T
       FDIP2=C*Q1
       RETURN
      END
C***********************************************************************
      REAL*8 FUNCTION FDIPA(EK1,L1,EK2,L2)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C N.R. BADNELL
C ASYMPTOTIC EXPRESSION FOR I(KAPPA1,L1,KAPPA2,L2,1) BASED ON A40,1 OF BHT
C
      IF(EK1*EK2.GT.1.0D-50)THEN
      X1=1.0D0/SQRT(EK1)
      X2=1.0D0/SQRT(EK2)
      XP=ABS(X1-X2)
      PI=ACOS(-1.0D0)
      XP=EXP(0.5D0*PI*XP)
      IF(EK1-EK2)1,1,2
   1  E=EK1/EK2
      IF(L1-L2)3,3,4
   3  L=L1
      GO TO 7
   4  L=L2
      GO TO 8
   2  E=EK2/EK1
      IF(L1-L2)5,5,6
   5  L=L1
      GO TO 8
   6  L=L2
      GO TO 7
C A40
   7  TL=L
      T0=1.0D0-E
      IF(TL*T0.LT.E)GO TO 9
      T=PI*TL
      EE=SQRT(E)
      F0=SQRT(T*T0*EE)*EE**L
      TL=L+L+1
      FDIPA=F0*XP/TL
      RETURN
C A41
   8  T0=1.0D0-E
      TL=L
      IF(TL*T0.LT.E)GO TO 9
      T0=1.0D0/T0
      T=TL*PI
      EE=SQRT(E)
      F0=SQRT(T*T0*EE)*EE**(L+1)
      TL=L+L+1
      TL2=L+L+3
      FDIPA=F0*XP/(TL*TL2)
      RETURN
      ENDIF
   9  FDIPA=0.0D0
      RETURN
      END
C
C***************************************************************
C
      REAL*8 FUNCTION FKHI(E,L,AC)
C
C  CALCULATES REAL*4 PART OF PSI(L+1+I*GAM) - LN(GAM)
C  WHERE E = 1/(GAM**2).
C  THIS IS REQUIRED FOR CALCULATION OF SCRIPT G.
C
      IMPLICIT REAL*8 (A-H,O-Y)
      IMPLICIT COMPLEX*16 (Z)
C
      PARAMETER (TZERO=0.0)
      PARAMETER (ONE=1.0)
      PARAMETER (TWO=2.0)
      PARAMETER (P0=ONE/252.)
      PARAMETER (P1=1.05)
      PARAMETER (P2=2.1)
C
      FKHI=TZERO
      IF(E.EQ.0)RETURN
C
      AC1=(20.*AC)**.333
C
      IF(E.GT.AC1)GOTO 100
C
      C=TZERO
      IF(L.EQ.0)GOTO 20
      A1=ONE
      A2=-E
      A3=E+E
      DO 10 I=1,L
      A2=A2+A3
      A1=A1+A2
      C=C+DBLE(I)/A1
   10 CONTINUE
   20 FKHI=E*((((P1*E+ONE)*E+P2)*E+21)*P0+C)
      RETURN
C
  100 AC1=ONE/SQRT(AC1)
      FL=DBLE(L+1)
      IF(FL.GT.AC1)GOTO 300
C
      N=AC1
      FL=N+1
      L1=L+1
      DO 210 I=L1,N
      FI=I
      FKHI=FKHI+FI/(ONE+E*FI*FI)
  210 CONTINUE
      FKHI=-FKHI*E
C
  300 X1=FL*E
      X=ONE+X1*FL
      ZE=DCMPLX(FL,ONE/SQRT(E))
      ZE=-ONE/(ZE*ZE)
      FKHI=FKHI+(LOG(X)-(X1/X))/TWO+DBLE((((P1*ZE+ONE)*ZE
     C +P2)*ZE+21)*ZE)*P0
C
      RETURN
      END
C***********************************************************************
       REAL*8 FUNCTION FMON1(EK1,EK2,L)
C
       IMPLICIT REAL*8 (A-H,O-Z)
C
       IF(EK1+EK2-1.0D-40)28,28,29
   28  FMON1=1.0D+50
       RETURN
   29  CONTINUE
       VMAX=200.0D0
       X1=SQRT(EK1)
       X2=SQRT(EK2)
       X3=X1+X2
       X4=X3*X3
       X5=X1*X2
       X6=X2-X1
       X7=4.0D0/X4
       PI=3.141592653589793D0
       IF(EK1-EK2)1,1,2
    1  ETA=1.0D0/X2
       GO TO 3
    2  ETA=1.0D0/X1
    3  G=0.5D0*PI*EXP(-PI*ETA)
       IF(G.EQ.0.0D0)GO TO 20           !NRB
       A1=1.0D0
       A2=1.0D0
       MG=0
       MA1=0
       MA2=0
       M=-1
    4  M=M+1
       EM=M
       T=EM+EM+1.0D0
       G=G*X7/(T*(T+1.0D0))
       EMM=EM*EM
       A1=A1*(1.0D0+EMM*EK1)
       A2=A2*(1.0D0+EMM*EK2)
   30  IF(G-0.015625D0) 31,32,32
   31  G=64.0D0*G
       MG=MG-1
       GO TO 30
   32  IF(G-64.0D0) 34,34,33
   33  G=0.015625D0*G
       MG=MG+1
       GO TO 32
   34  IF(A1-64.0D0) 36,36,35
   35  A1=0.015625D0*A1
       MA1=MA1+1
       GO TO 34
   36  IF(A2-64.0D0) 38,38,37
   37  A2=0.015625D0*A2
       MA2=MA2+1
       GO TO 36
   38  CONTINUE
       IF(M-L)4,5,5
    5  G=G*(T+1.0D0)
       IF(X1-300.0D0)7,6,6
    6  B=PI/X1
       A1=1.5D0*A1/(B*(3.0D0-B*(3.0D0-B*(2.0D0-B))))
       GO TO 9
    7  IF(X1-0.2D0)9,9,8
    8  B=-PI/X1
       A1=A1/(1.0D0-EXP(B+B))
    9  IF(X2-300.0D0)11,10,10
   10  B=PI/X2
       A2=1.5D0*A2/(B*(3.0D0-B*(3.0D0-B*(2.0D0-B))))
       GO TO 13
   11  IF(X2-0.2D0)13,13,12
   12  B=-PI/X2
       A2=A2/(1.0D0-EXP(B+B))
   13  G=G*SQRT(A1*A2)*(8.0D0)**(MG+MG+MA1+MA2)
       S0=1.0D0
       S1=0.0D0
       U=L
       V=0.0D0
       W=U+U+1.0D0
       T0=1.0D0
       T1=0.0D0
   14  U=U+1.0D0
       V=V+1.0D0
       W=W+1.0D0
       IF(V-VMAX)21,21,20
   20  FMON1=0.0D0
       RETURN
   21  CONTINUE
       U0=U*U*X5+1.0D0
       U1=U*X6
       T=T0*U0-T1*U1
       T1=T0*U1+T1*U0
       T0=T
       T=X7/(V*W)
       T0=T*T0
       T1=T*T1
       S0=S0+T0
       S1=S1+T1
       S=S0*S0+S1*S1
       T=T0*T0+T1*T1
       SM=1.0D0/S
       TM=1.0D0/T
       IF(SM*TM.EQ.0.0D0)GO TO 20      !NRB
       IF(S-1.0D+24*T)14,15,15
   15  FMON1=G*SQRT(S)
       IV=V
       RETURN
      END
C***************************************************************

C***********************************************************************
       REAL*8 FUNCTION ARGAM(L,A)
C
       IMPLICIT REAL*8 (A-H,O-Z)
C
C CALCULATES ARGGAMMA(L+1+I*A)
C WHERE L IS AN INTEGER NOT LESS THAN ZERO
C
       B=ABS(A)
       B=250.0D0*B**0.25D0-A*A
       J0=L+1
       C=J0
       D=C*C
       Z=0.0D0
       IF(D -B)1,6,6
    1  B=SQRT (B)
       J1=B
       DO 5 J=J0,J1
       D=J
       D=A/D
       G1=ABS(D)
       IF(G1-0.1D0)2,3,3
    2  G1=D*D
       G2=-35.0D0*G1+45.0D0
       G2=-G1*G2+63.0D0
       G2=-G1*G2+105.0D0
       G1=D -D*G1*G2/315.0D0
       GO TO 4
    3  G1=ATAN (D)
    4  Z=Z+G1
    5  CONTINUE
       J0=J1+1
    6  D=J0
       G0=D*D
       U=A*A
       G1=1.0D0/(G0+U)
       G2=G1*G1
       G3=10.0D0*G0*G0-20.0D0*G0*U+2.0D0*U*U
       G3=G3*G2-21.0D0*G0+7.0D0*U
       G3=G3*G2+210.0D0
       G1=A*G3*G1/2520.0D0
        ARGAM=-Z+0.5D0*A*LOG(G0+U)+(D -0.5D0)*ATAN(A/D)-A-G1
       RETURN
      END
C***************************************************************
