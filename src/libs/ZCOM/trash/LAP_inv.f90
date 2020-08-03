!======================================================================
      Subroutine LAP_INV(n,m,A,info)
!======================================================================
!     Call LAPACK procedures DGETRF and DGETRI to obtain the inverse
!     matrix A^-1 for symmetric matrix A
!---------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(4), INTENT(in) :: n,m
      REAL(8), DIMENSION(m,m) :: A
      Integer, intent(out), optional :: info

      REAL(8), ALLOCATABLE, DIMENSION(:) :: work
      INTEGER(4), ALLOCATABLE, DIMENSION(:) :: ipiv
      INTEGER(4) :: lwork,ierr

      lwork = 3*n
      Allocate(IPIV(n), WORK(lwork), stat=ierr)
      if(ierr.ne.0) Stop 'LAP_INV: problems with allocation '

      Call DGETRF(n,n,A,m,IPIV,info)

      if(info.ne.0) write(*,*) 'LAP_INV: DGETRF give INFO =',INFO

      Call DGETRI(n,A,m,IPIV,WORK,LWORK,INFO)

      if(info.ne.0) write(*,*) 'LAP_INV: DGETRI give INFO =',INFO

      Deallocate(IPIV, WORK)

      End Subroutine LAP_INV


!======================================================================
      Subroutine LAP_INVS(UPLO,n,m,A,info)
!======================================================================
!
!     Call LAPACK procedures DGETRF and DGETRI to obtain the inverse
!     matrix A^-1 for symmetric matrix 
!
!---------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(4), INTENT(in) :: n,m
      Character(1), Intent(in) :: UPLO
      REAL(8), DIMENSION(m,m) :: A
      Integer, intent(out), optional :: info

      REAL(8), ALLOCATABLE, DIMENSION(:) :: work
      INTEGER(4), ALLOCATABLE, DIMENSION(:) :: ipiv
      INTEGER(4) :: lwork,ierr

      lwork = 3*n
      Allocate(IPIV(n), WORK(lwork), stat=ierr)
      if(ierr.ne.0) Stop 'LAP_INV: problems with allocation '

      Call DSYTRF(UPLO,n,A,m,IPIV,WORK,LWORK,INFO)

      if(info.ne.0) write(*,*) 'LAP_INV: DGETRF give INFO =',INFO
      if(info.ne.0) Return

      Call DSYTRI(UPLO,n,A,m,IPIV,WORK,INFO)

      if(info.ne.0) write(*,*) 'LAP_INV: DGETRI give INFO =',INFO

      Deallocate(IPIV, WORK)

      End Subroutine LAP_INVS
