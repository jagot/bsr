!======================================================================
      Real(8) Function GEN_INT(N,IA)
!======================================================================
!     estimates the gen = N * Product [ P(j)^IA(j) ]
!----------------------------------------------------------------------
      Implicit none
      Integer(2) :: N
      Integer(1) :: IA(11)
      Real(8) :: C,P(11)
      Data P/2.d0,3.d0,5.d0,7.d0,11.d0,13.d0,17.d0,19.d0,23.d0,29.d0,31.d0/ 
      Integer :: i,j,m

      C=0.d0
      if(N.ne.0) then
       C=1.d0
       Do J=1,11
        M=abs(IA(J))
        if(IA(J).gt.0) then
         Do i=1,M; C=C*P(J); End do
        end if
        if(IA(J).lt.0) then
         Do i=1,M; C=C/P(J); End do
        end if
       End do
      end if

      GEN_INT = DFLOAT(N)*dsqrt(C)

      End Function GEN_INT
