!---------------------------------------------------------------------
      Real(8) FUNCTION GRAD3(I,J)
!---------------------------------------------------------------------
!
!     computes the dipole integral between i-th and j-th orbitals
!     in the V-form:  INT{ P_j(r) [r d/dr + f(l)] P_i(r); dr}
!
!---------------------------------------------------------------------

      Use RADIAL, L => lro

      Implicit real(8) (A-H,O-Z)
      Integer(4), Intent(in) :: I,J

      DIMENSION Q(NR)

      LI=L(I); LJ=L(J); IL = IABS(LI - LJ)
      IF (IL .NE. 0 .AND. IL .NE. 2) Stop 'GRAD2: LI - LJ <> 0,2'

      IP = nrf+1;  Call P_derive(J,IP)

      G = QUADR(I,IP,1); U = QUADR(I,J,0)

      idelta=LJ-LI

      Select case(idelta)
       Case(-2); GRAD3=G-(LI-D2)*U
       Case( 0); GRAD3=G+(D1+D5)*U
       Case(+2); GRAD3=G+(LI+D3)*U
      End Select

      END FUNCTION GRAD3

