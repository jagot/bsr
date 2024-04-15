!======================================================================
      Subroutine LAP_DSYEV(job,UPLO,n,m,A,eval,info)
!======================================================================
!
!     Call LAPACK procedure DSYEV to obtain the eigenvalues (eval) and
!     eigenvectors (A) for real symmetric matrix A(n,n)
!
!     job = 'V' or 'N' - compute or not the eigenvectors
!
!     used upper UPLO='U' triangle matrix, or lower UPLO='L'
!---------------------------------------------------------------------

      IMPLICIT NONE

      Character(1), INTENT(in) :: job,UPLO
      INTEGER(4), INTENT(in) :: n,m
      Integer, intent(out), optional :: info
      REAL(8), DIMENSION(m,m) :: A
      REAL(8), DIMENSION(m) :: eval

      REAL(8), ALLOCATABLE, DIMENSION(:) :: work
      INTEGER(4) :: lwork

      lwork = 3*n-1; Allocate(work(lwork))

      Call DSYEV(job,UPLO,n,A,m,eval,WORK,LWORK,INFO )

      if(INFO.ne.0) write(*,*) 'LAP_DSYEV: DSYEV(lapack) gives INFO = ',INFO

      Deallocate(work)

      End Subroutine LAP_DSYEV
