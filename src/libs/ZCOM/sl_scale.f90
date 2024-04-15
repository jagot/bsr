!======================================================================
      Subroutine SL_scale(N,ro1,h,A,B,R,RT,acc)
!======================================================================
!     Semi-logarithmic scale:   ro = b*ln(r) + a*r
!     ro is defined at N points:  ro1 + (i-1)*h,  i=1,N
!     RT - initial estimation
!     acc - accuracy for Newton method
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in)  :: N
      Real(8), intent(in)  :: ro1,h,A,B,acc
      Real(8), intent(out) :: R(N)
      Integer :: i
      Real(8) :: RT,RO,S

      DO I=1,N
       RO=ro1+(I-1)*h
       DO
        S=(A*RT+B*LOG(RT)-RO)/(A*RT+B)
        RT=RT*(1.D0-S)
        IF(ABS(S).lt.acc) Exit
       END DO
       R(I)=RT
      END DO

      End Subroutine SL_scale
