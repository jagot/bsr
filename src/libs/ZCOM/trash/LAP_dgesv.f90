!======================================================================
      Subroutine LAP_DGESV(m,n,k,A,B,info)
!======================================================================
!     Call LAPACK procedure DGESV to solve the system of algebraic
!     equations A x = B, results in B.
!---------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(4), INTENT(in) :: n,m,k
      REAL(8), DIMENSION(m,*) :: A,B
      Integer, intent(out), optional :: info

      Integer(4), ALLOCATABLE, DIMENSION(:) :: IPIV

      Allocate(IPIV(m))

      Call DGESV( n,k, A, m, IPIV, B, m, INFO )

      if(info.ne.0) write(*,*) ' DGESV(lapack) give INFO =',INFO

      Deallocate(IPIV)

      End Subroutine LAP_DGESV
