!======================================================================
      Subroutine LAP_DSYGV(job,UPLO,n,m,A,C,eval,info)
!======================================================================
!
!     Call LAPACK procedure DSYEV to obtain the eigenvalues (eval) and 
!     eigenvectors (A) for generalized problem A S = E C S
!
!     job = 'N' - compute eigenvalues only
!           'V' - compute eigenvalues and eigenvectors
!
!     used upper UPLO='U' triangle matrix, or lower UPLO='L'
!---------------------------------------------------------------------

      IMPLICIT NONE

      Character(1), INTENT(in) :: job,UPLO
      INTEGER(4), INTENT(in) :: n,m
      Integer, intent(out), optional :: info
      REAL(8), DIMENSION(m,*) :: A,C
      REAL(8), DIMENSION(*) :: eval

      INTEGER(4) :: lwork
      REAL(8), ALLOCATABLE, DIMENSION(:) :: work

      lwork = 3*n;  Allocate(work(lwork))

      Call DSYGV(1,job,UPLO,n,A,m,C,m,eval,WORK,LWORK,INFO)

      if(INFO.ne.0) write(*,*) ' DSYGV(lapack) gives INFO = ',INFO

      Deallocate(work)

      End Subroutine LAP_DSYGV
