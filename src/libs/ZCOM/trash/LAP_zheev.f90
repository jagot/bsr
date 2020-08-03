!======================================================================
      Subroutine LAP_ZHEEV(job,UPLO,n,m,A,eval,info)
!======================================================================
!     Call LAPACK procedure ZHEEV to obtain the eigenvalues (eval) and 
!     eigenvectors (A) for complex Hermitian matrix A(n,n)
!
!     job = 'V' or 'N' - compute or not the eigenvectors
!
!     used upper UPLO='U' triangle matrix, or lower UPLO='L'
!---------------------------------------------------------------------
      IMPLICIT NONE
      Character(1), INTENT(in) :: job,UPLO
      INTEGER(4), INTENT(in) :: n,m
      Integer, intent(out), optional :: info
      COMPLEX(8), DIMENSION(m,m) :: A
      REAL(8), DIMENSION(m) :: eval

      REAL(8), ALLOCATABLE :: rwork(:)
      COMPLEX(8), ALLOCATABLE :: work(:)
      INTEGER(4) :: lwork,lrwork

      lwork = 2*n; lrwork= 3*n
      Allocate (WORK(lwork), RWORK(lrwork))

      Call ZHEEV(job,UPLO,n,A,m,eval,WORK,lwork,RWORK,INFO)

      if(INFO.ne.0) write(*,*) ' ZHEEV(lapack) gives INFO = ',INFO

      Deallocate(work,rwork)

      End Subroutine LAP_ZHEEV
