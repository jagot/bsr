!=====================================================================
      Real(8) Function XLAGR(k,n,X,Y,xx)
!=====================================================================
!     Interpolation with Lagrang formula for 1 point (xx)
!     based on k points in function Y(x).
!     XLAGR = 0 if xx is outside the interval X(:)
!
!     Call: INTERV, RLAGR
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k,n
      Real(8), intent(in) :: X(n),Y(n)
      Real(8), intent(in) :: XX
      Integer :: left, mflag
      Real(8), external :: RLAGR 

      if(k.gt.n) Stop ' XLAGR: K > N'

      Call INTERV ( x, n, xx, left, mflag )

      if(mflag.ne.0) then
       XLAGR=0.0
      end if

!      left = max(1,left-k/2)
      if(left+k-1.gt.n) left = n-k+1

      XLAGR = RLAGR(k, X(left), Y(left), xx)

      End function XLAGR


!===================================================================== 
      Real(8) Function RLAGR (K,X,Y,XX) 
!===================================================================== 
!     Interpolation with Lagrang formula (on K points) for 1 point (XX) 
!--------------------------------------------------------------------- 
      Implicit none 
      Integer, intent(in) :: k
      Real(8), intent(in) :: X(k),Y(k),XX
      Real(8) :: YY
      Integer :: I,J
 
      rlagr = 0.d0 
 
      DO J=1,K 
        YY = Y(j) 
        DO I=1,K 
          IF(I.NE.J) YY=YY*(XX-X(i))/(X(j)-X(i)) 
        END DO 
        rlagr = rlagr + yy 
      END DO 
 
      End Function RLAGR


!====================================================================== 
      Subroutine LAGRN(K,N,N1,R,R1,F,F1) 
!====================================================================== 
!     Lagrang interpolation  F(R) --> F1(R1) 
!
!     Call: XLAGR 
!---------------------------------------------------------------------- 
      Implicit none
      Integer, intent(in ) :: K,N,N1
      Real(8), intent(in ) :: R(N),F(N)
      Real(8), intent(out) :: R1(N1),F1(N1)
      Integer :: I1
      Real(8), external :: XLAGR

      Do I1=1,N1 
       F1(i1) = XLAGR(k,n,R,F,R1(i1)) 
      End do 
 
      End Subroutine LAGRN



