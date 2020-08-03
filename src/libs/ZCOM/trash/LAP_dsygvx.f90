!======================================================================
      Subroutine LAP_DSYGVX(job,UPLO,n,m,A,C,eval,k,info)
!======================================================================
!     Call LAPACK procedure DSYGVX to obtain first k eigenvalues 
!     (eval) and eigenvectors (A) for generalized problem  A S = E C S
!
!     job = 'N' - compute eigenvalues only
!           'V' - compute eigenvalues and eigenvectors
!
!     UPLO = 'U' - used upper triangle matrix
!            'L' - used lower triangle matrix
!---------------------------------------------------------------------
      IMPLICIT NONE
      Character(1), INTENT(in) :: job,UPLO
      INTEGER(4), INTENT(in) :: n,m,k
      Integer, intent(out), optional :: info
      REAL(8), DIMENSION(m,*) :: A,C
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

      Call DSYGVX(1,job,'I',UPLO,n,A,m,C,m,VL,VU,1,k,ABSTOL,kk,eval, &
                  Z,n, WORK, LWORK, IWORK, IFAIL, INFO)
      
      if(kk.ne.k) then
       write(*,*) ' DSYGVX(lapack) provides',kk,' eigenvectors'
       write(*,*) ' when we ordered', k
      end if
      if(INFO.ne.0) write(*,*) 'DSYGVX(lapack) gives INFO = ',INFO

      A(1:n,1:k) = Z(1:n,1:k)

      Deallocate(work,ifail,Z)

      End Subroutine LAP_DSYGVX
