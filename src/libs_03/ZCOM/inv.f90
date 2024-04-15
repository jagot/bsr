!======================================================================
      Subroutine INV (N,mdim,A)
!======================================================================
!     Inverse matrix A(n,n) (from IMPACT code)
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Integer, intent(in) :: N, mdim
      Real(8), intent(inout) :: A(mdim,mdim)
      Integer  :: IP(N), IND1(N), IND2(N)

      if(n.eq.1) then
        A(1,1)=1./A(1,1)
        Return
      elseif(n.le.0) then
        Stop ' INV: n <= 0 '
      elseif(n.gt.mdim) then
        Stop ' INV: n > mdim '
      end if

      IP = 0

      DO I = 1,N
       AMAX = 0.0

       DO J = 1,N
        IF(IP(J).eq.1) Cycle
        DO K = 1,N
         IF(IP(K).eq.1) Cycle
         IF(IP(K).gt.1) Return
         IF(ABS(AMAX).ge.ABS(A(J,K))) Cycle
         IR = J
         IC = K
         AMAX = A(J,K)
        END DO
       END DO

       IP(IC) = IP(IC) + 1

       if(IR.ne.IC) then
        Do L = 1,N
         SWAP = A(IR,L)
         A(IR,L) = A(IC,L)
         A(IC,L) = SWAP
        End do
       end if

       IND1(I) = IR
       IND2(I) = IC

       PIVI = 1./A(IC,IC)
       A(IC,IC) = 1.
       Do L = 1,N
        A(IC,L) = A(IC,L) * PIVI
       End do

       DO L1 = 1,N
        IF(L1.eq.IC) Cycle
        T = A(L1,IC)
        A(L1,IC) = 0.
        Do L = 1,N
         A(L1,L) = A(L1,L) - A(IC,L)*T
        End do
       End do

      End do  !  main loop over I

      L = N + 1
      Do I=1,N
       L = L - 1
       if(IND1(L).eq.IND2(L)) Cycle
       JR = IND1(L)
       JC = IND2(L)
       Do K=1,N
        SWAP=A(K,JR)
        A(K,JR)=A(K,JC)
        A(K,JC)=SWAP
       End do
      End do

      End Subroutine INV



