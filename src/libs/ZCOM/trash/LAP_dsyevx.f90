!======================================================================
      Subroutine LAP_DSYEVX(job,UPLO,n,m,A,eval,k,info)
!======================================================================
!
!     Call LAPACK procedure DSYEVX to obtain the selected eigenvalues 
!     (eval) and eigenvectors (A) for problem A S = E S
!
!     used upper UPLO='U' triangle matrix, or lower UPLO='L'
!---------------------------------------------------------------------

      IMPLICIT NONE

      Character(1), INTENT(in) :: job,UPLO
      INTEGER(4), INTENT(in) :: n,m,k
      Integer, intent(out), optional :: info
      REAL(8), DIMENSION(m,*) :: A
      REAL(8), DIMENSION(*) :: eval

      REAL(8), ALLOCATABLE, DIMENSION(:) :: work
      REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Z
      INTEGER(4), ALLOCATABLE, DIMENSION(:) :: iwork,ifail

      INTEGER(4) :: kk, lwork
      REAL(8) :: ABSTOL, VL =0.d0, VU = 0.d0
      REAL(8), EXTERNAL :: DLAMCH

      ABSTOL = 4*DLAMCH('S')
      lwork = 10*n
      Allocate(work(lwork),iwork(5*n),ifail(k),Z(n,k))

      Call DSYEVX(job,'I',UPLO,n,A,m,VL,VU,1,k,ABSTOL,kk,eval, &
                  Z,n, WORK, LWORK, IWORK, IFAIL, INFO)
      
      if(kk.ne.k) then
       write(*,*) ' DSYEVX(lapack) provides',kk,' eigenvectors'
       write(*,*) ' when we ordered', k,'  ones'
       Stop ' Stop in LAP_DSYEVX'
      end if

      if(INFO.ne.0) write(*,*) 'DSYEVX(lapack) gives INFO = ',INFO

      A(1:n,1:k) = Z(1:n,1:k)

      Deallocate(work,ifail,Z)

      End Subroutine LAP_DSYEVX
