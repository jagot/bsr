!------------------------------------------------------------------
      SUBROUTINE DYK(k,kp)
!------------------------------------------------------------------
!
!     Computes the  solution of the differential equation
!
!           (d - (k+kp)/r) yk(r)= -(2k+kp)(1/r) zk(r)
!
!     with yk(no) = zk(no), yk(no-1) = zk(no-1)
!
!     kp = 1  -->  RK, TK integrals
!     kp = 3  -->  NK, VK integrals
!
!     Uses the inward integration with finate-difference formula.
!
!     INPUT  zk(r)  -   ZK array in RADIAL
!     OUTPUT yk(r)  -   YK array in RADIAL
!------------------------------------------------------------------

      USE RADIAL

      Implicit real(8) (A-H,O-Z)

      Integer(4), Intent(in) :: k,kp

      B  = EH**k
      A  = B*EH**kp
      AA = A*A
      A  = 4*A
      HH = (k+k+kp)*H/D3
      F2 = ZK(NR-2)*B
      F1 = F2*B
      Do M = NR-2,1,-1
       F3 = ZK(M)
       YK(M) = YK(M+2)*AA + HH*(F3 + A*F2 + AA*F1)
       F1 = F2
       F2 = F3
      End do

      END SUBROUTINE DYK
