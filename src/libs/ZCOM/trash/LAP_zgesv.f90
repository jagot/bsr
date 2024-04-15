!======================================================================
      Subroutine LAP_ZGESV(n,m,A,B,info)
!======================================================================
!
!     Call LAPACK procedure ZGESV to solve the system of complex
!     algebraic  equations A x = B, results in B.
!
!---------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(4), INTENT(in) :: n,m
      COMPLEX(8), DIMENSION(m,m) :: A,B
      Integer, intent(out), optional :: info

      Integer(4), ALLOCATABLE, DIMENSION(:) :: IPIV

      Allocate(IPIV(m))

      Call ZGESV( m,n, A, m, IPIV, B, m, INFO )

      if(info.ne.0) write(*,*) ' ZGESV(lapack) give INFO =',INFO

      Deallocate(IPIV)

      End Subroutine LAP_ZGESV
