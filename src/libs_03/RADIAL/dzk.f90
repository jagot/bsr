!------------------------------------------------------------------
      SUBROUTINE dzk (I,J,K)
!------------------------------------------------------------------
!
!     Solves the following differential equation
!
!        (d - k/r) zk(r) = P_i(r) * P_j(r),
!
!     with boundary conditionas zk(0) = 0
!
!     zk(1),zk(2),zk(3) are approximated by expansions in r
!
!     Result:  array ZK (module RADIAL)
!              zk(r) ~ r^m (1 + a*(rz) + b*(rz)^2 + ...)
!
!------------------------------------------------------------------

      USE RADIAL

      IMPLICIT NONE
      Integer(4), Intent(in) :: i,j,k
      Integer(4) :: m, m1,m2
      Real(8) :: A,B,C, A1,A2,A3,A4,A5, F1,F2,F3,F4,F5

! ... estimation of ZK in first points ...

      m = mexp(i) + mexp(j) + k
      C = m + 1
      A = aexp(i) + aexp(j)
      B = bexp(i) + bexp(j) + aexp(i)*aexp(j)
      B = A*A/(C+D1) - D2*B/(C+D2)
      A = -A/(C+D1)

      F1 = RR(1)*P(1,I)*P(1,J)
      F2 = RR(2)*P(2,I)*P(2,J)
      F3 = RR(3)*P(3,I)*P(3,J)
      F4 = RR(4)*P(4,I)*P(4,J)

      ZK(1) = F1*(D1 + Z*R(1)*(A + Z*R(1)*B) )/C
      ZK(2) = F2*(D1 + Z*R(2)*(A + Z*R(2)*B) )/C
      ZK(3) = F3*(D1 + Z*R(3)*(A + Z*R(3)*B) )/C

! ... power expantion for Zk ...

      mzk = mexp(i) + mexp(j) + 1
      azk = (aexp(i) + aexp(j))*C/(C+D1)
      bzk = (bexp(i) + bexp(j) + aexp(i)*aexp(j))*C/(C+D2)

! ... integration of equation ...

      A   =  EH**K
      B   =  H/90.D0
      C   =  A*A
      A1  = -C*A*B
      A2  =  34.D0*C*B
      A3  =  114.D0*A*B
      A4  =  34.D0*B
      A5  = -B/A

      DO M = 5,NR
       F5 = RR(M)*P(M,I)*P(M,J)
       ZK(M-1) = ZK(M-3)*C + F1*A1 + F2*A2 + F3*A3 + F4*A4 + F5*A5
       F1 = F2; F2 = F3; F3 = F4; F4 = F5
      End do
      ZK(NR) = ZK(NR-1) * A

!     FOR Z0(I,I) SET THE LIMIT TO 1 AND REMOVE OSCILLATIONS
!     INTRODUCED BY THE USE OF SIMPSON'S RULE

      if (IABS(I-J) + IABS(K).eq.0) then
       M1 = (NR/2)*2 - 1
       M2 = M1 - 1
       F1 = D1 - ZK(M1)
       F2 = D1 - ZK(M2)
       Do M = 1,M1,2
        ZK(M  ) = ZK(M  ) + F1
        ZK(M+1) = ZK(M+1) + F2
       End do
       ZK(NR  ) = D1
       ZK(NR-1) = D1
      end if

      END SUBROUTINE dzk
