!---------------------------------------------------------------------
      Real(8) Function DET(N,A)
!---------------------------------------------------------------------
!     determinant of array A(N,N)   (Gauss method)
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: N
      Real(8), intent(inout) :: A(N,N)
      Integer :: I,J,K
      Real(8) :: MAX, T

      DET = 1.d0

      DO K=1,N

       MAX=0.d0
       DO I=K,N
        T=A(I,K)
        if(ABS(T).gt.ABS(MAX)) then
         MAX=T; J=I
        end if
       END DO

       IF(MAX.EQ.0.d0) THEN; DET=0.d0; Return; END IF

       IF(J.NE.K) THEN
        DET = -DET
        DO I=K,N
         T=A(J,I); A(J,I)=A(K,I); A(K,I)=T
        END DO
       END IF

       IF(K+1.LE.N) THEN
        DO I=K+1,N
         T=A(I,K)/MAX
         DO J=K+1,N
          A(I,J)=A(I,J)-T*A(K,J)
         END DO
        END DO
       END IF

       DET=DET*A(K,K)

      END DO

      End Function DET


!=====================================================================
      Real(8) Function DETA(N,M,A)
!=====================================================================
!     determinant of matrix A(N,N)   (Gauss method)
!=====================================================================
      Implicit none
      Integer, intent(in) :: N,M
      Real(8), intent(inout) :: A(M,*)
      Integer :: I,J,K
      Real(8) :: MAX, T

      DETA = 1.d0

      DO K=1,N

       MAX=0.d0
       DO I=K,N
        T=A(I,K)
        if(ABS(T).gt.ABS(MAX)) then
         MAX=T; J=I
        end if
       END DO

       IF(MAX.EQ.0.d0) THEN; DETA=0.d0; Return; END IF

       IF(J.NE.K) THEN
        DETA = -DETA
        DO I=K,N
         T=A(J,I); A(J,I)=A(K,I); A(K,I)=T
        END DO
       END IF

       IF(K+1.LE.N) THEN
        DO I=K+1,N
         T=A(I,K)/MAX
         DO J=K+1,N
          A(I,J)=A(I,J)-T*A(K,J)
         END DO
   	    END DO
       END IF

       DETA=DETA*A(K,K)

      END DO

      End Function DETA
