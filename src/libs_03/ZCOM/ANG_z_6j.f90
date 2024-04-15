!======================================================================
      Real(8) Function Z_6J (j1,j2,j3,j4,j5,j6)
!======================================================================
!     determination of 6j-symbols without direct using of factorials
!     accoding to formula:
!
!     6j{j1,j2,j3,j4,j5,j6) = {j1,j2,j3}*{j1,j5,j6}*{j4,j2,j3}*{j4,j5,j3}*
!                                SUM(z) {   (-1)^z * (z+1)!   /
!          [ (z-j1-j2-j3)! * (z-j1-j5-j6)! * (z-j4-j2-j3)! *(z-j4-j5-j3)! *
!              (j1+j2+j4+j5-z)! * (j1+j3+j4+j6-z)! * (j2+j3+j5+j6-z)! ]
!
!     where {a,b,c}=sqrt[ (a+b-c)! * (a-b+c)! * (b-a+c)! / (a+b+c+1)! ]
!
!     If we introduce the auxiliary values L(i)
!     (see below the text of program) then
!
!     6j = sqrt{ Pr(j=5,7,i=1,4) (L(j)-L(i))! / Pr(i=1,4) (L(i)+1)! }
!                Sum(z) { (-1)^z * (z+1)! /
!          [ Pr(i=1,4) (z-L(i))!  * Pr(j=5,7) (L(j)-z) ] }
!
!     (the momenta are used in (2J+1)-representation)
!
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j1,j2,j3,j4,j5,j6
      Integer :: I, IZ_min, IZ_max, K, KK, M, IZ, I_max
      Integer :: L(7),J(23)
      Real(8) :: X, R, C

      Z_6J = 0.0

      L(1)=j1+j2+j3-3                    ! auxiliary values
      L(2)=j1+j5+j6-3
      L(3)=j4+j2+j6-3
      L(4)=j4+j5+j3-3
      L(5)=j1+j2+j4+j5-4
      L(6)=j1+j3+j4+j6-4
      L(7)=j2+j3+j5+j6-4
      DO I=1,7
      IF(mod(L(I),2).eq.1) Return
       L(I)=L(I)/2
      END DO

      IZ_min=MAX0(L(1),L(2),L(3),L(4))   ! limits of the sum
      IZ_max=MIN0(L(5),L(6),L(7))
      IF(IZ_max.LT.IZ_min) Return

      Do I=1,4
       J(I)=L(I)+1
       Do K=5,7
        M=L(K)-L(I)
        IF(M.LT.0) Return                ! check of triangle rule
        J(4+3*I+K)=M
       End do
      End do
                                         ! initial factorial parameters
                                         ! in the sum
      Do I=5,8;  J(I)=IZ_min-L(I-4); End do
      Do I=9,11; J(I)=L(I-4)-IZ_min; End do

      C=0.0
      DO IZ=IZ_min,IZ_max           ! summation
       I_max=IZ+1

!      this limit for max. factorial follows from symmetry propeties:
!      let's a(i)=L(i) for i=1,4;  b(i)=L(j),j=5,7;
!      then  b(j)-a(i) <= a(k)  <= max[a(k)] = IZ_min;
!      also  a(j) <= b(j), then a(i)-a(j) <= b(i)-a(j) <= IZ_min;
!      and last (a(i)+1) < max(a(i))+1 < IZ_min+1;

       X=1.0
       DO I=2,I_max            ! estimation of one term in sum
        K=2                    ! K - the power of integer I in the term
        DO M=12,23; IF(J(M).GE.I) K=K+1;  End do
        DO M=1,4;   IF(J(M).GE.I) K=K-1;  End do
        DO M=5,11;  IF(J(M).GE.I) K=K-2;  End do
        IF(K.EQ.0) Cycle

        R=DBLE(I)               ! X = X * I ** K/2
        KK=IABS(K)/2
        IF(KK.GT.0) THEN
         DO M=1,KK
          IF(K.GT.0) X=X*R
          IF(K.LT.0) X=X/R
         END DO
        END IF
        IF(mod(K,2).EQ.+1) X=X*SQRT(R)
        IF(mod(K,2).EQ.-1) X=X/SQRT(R)
       End do

       IF(mod(IZ,2).eq.1) X=-X
       C=C+X
                                    ! new factorial parameters in sum
       Do I=5,8;  J(I)=J(I)+1; End do
       DO I=9,11; J(I)=J(I)-1; End do

      End do                        ! end of summation

      Z_6J=C

      End Function Z_6J


!====================================================================
      Real(8) Function Z_6jj(j1,j2,j3,j4,j5,j6)
!====================================================================

      Implicit none
      Integer, intent(in) :: j1,j2,j3,j4,j5,j6
      Real(8), external :: Z_6j

      Z_6jj = Z_6j(j1+j1+1,j2+j2+1,j3+j3+1,j4+j4+1,j5+j5+1,j6+j6+1)

      End Function Z_6JJ

!====================================================================
      Real(8) Function Z_6j2(j1,j2,j3,j4,j5,j6)
!====================================================================

      Implicit none

      Integer, intent(in) :: j1,j2,j3,j4,j5,j6
      Real(8), external :: Z_6j

      Z_6j2 = Z_6j(j1+1,j2+1,j3+1,j4+1,j5+1,j6+1)

      End Function Z_6j2
