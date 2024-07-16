!----------------------------------------------------------------------
!  interface routines to call some frequently-used LAPACK routines
!
!  LAP_DGESV     -    A x = B,      A,B - real
!  LAP_ZGESV     -    A x = B,      A,B - complex
!  LAP_DSYEV     -    A x = E x,    A - real symmetric
!  LAP_DSYEVX    -    A x = E x,    for some eigenvalues only
!  LAP_ZHEEV     -    A x = E x,    A - complex Hermitian
!  LAP_DSYGV     -    A x = E C x,  A,C - real symmetric
!  LAP_DSYGVX    -    A x = E C x,  for some eigenvalues only
!  LAP_INV       -    A^-1,         A - real
!  LAP_INVS      -    A^-1,         A - real symmetric
!----------------------------------------------------------------------

!======================================================================
      Subroutine LAP_DGESV(m,n,k,A,B,info)
!======================================================================
!     Call LAPACK procedure DGESV to solve the system of algebraic
!     equations A x = B, results in B.
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,m,k
      Real(8) :: A(m,*),B(m,*)
      Integer, intent(out), optional :: info
      Integer :: IPIV(m)

      Call DGESV( n,k, A, m, IPIV, B, m, INFO )

      if(info.ne.0) write(*,*) ' DGESV(lapack) give INFO =',INFO

      End Subroutine LAP_DGESV


!======================================================================
      Subroutine LAP_ZGESV(n,m,A,B,info)
!======================================================================
!     Call LAPACK procedure ZGESV to solve the system of complex
!     algebraic  equations A x = B, results in B.
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,m
      COMPLEX(8) :: A(m,m),B(m,m)
      Integer, intent(out), optional :: info
      Integer :: IPIV(m)

      Call ZGESV( m,n, A, m, IPIV, B, m, INFO )

      if(info.ne.0) write(*,*) ' ZGESV(lapack) give INFO =',INFO

      End Subroutine LAP_ZGESV


!======================================================================
      Subroutine LAP_DSYEV(job,UPLO,n,m,A,eval,info)
!======================================================================
!     Call LAPACK procedure DSYEV to obtain the eigenvalues (eval) and
!     eigenvectors (A) for Real symmetric matrix A(n,n)
!     job = 'V' or 'N' - compute or not the eigenvectors
!     UPLO='U'  or 'L' -  used upper or lower triangle matrix
!---------------------------------------------------------------------
      Implicit none
      Character(1), intent(in) :: job,UPLO
      Integer, intent(in) :: n,m
      Integer, intent(out) :: info
      Real(8) :: A(m,m)
      Real(8) :: eval(m)
      Real(8) :: work(3*n-1)
      Integer :: lwork

      lwork = 3*n-1

      Call DSYEV(job,UPLO,n,A,m,eval,WORK,LWORK,INFO )

      if(INFO.ne.0) write(*,*) 'LAP_DSYEV: DSYEV(lapack) gives INFO = ',INFO

      End Subroutine LAP_DSYEV


!======================================================================
      Subroutine LAP_DSYEVX(job,UPLO,n,m,A,eval,k,info)
!======================================================================
!     Call LAPACK procedure DSYEVX to obtain the selected eigenvalues
!     (eval) and eigenvectors (A) for problem A S = E S
!     UPLO='U'  or 'L' -  used upper or lower triangle matrix
!---------------------------------------------------------------------
      Implicit none
      Character(1), intent(in) :: job,UPLO
      Integer, intent(in) :: n,m,k
      Integer, intent(out), optional :: info
      Real(8) :: A(m,*)
      Real(8) :: eval(*)

      Real(8), allocatable :: work(:), Z(:,:)
      Integer, allocatable :: iwork(:),ifail(:)
      Integer :: kk, lwork
      Real(8) :: ABSTOL, VL =0.d0, VU = 0.d0
      Real(8), external :: DLAMCH
      Character(256) :: error_msg

      if(m<n) then
         write(error_msg, '("m = ",i10," < n = ",i10, ", m >= n required")') m, n
         error stop trim(error_msg)
      end if

      if(k<1) then
         write(error_msg, '("k = ",i10," < 1, k >= 1 required")') k
         error stop trim(error_msg)
      end if

      ABSTOL = 4*DLAMCH('S')
      lwork = 10*n
      Allocate(work(lwork),iwork(5*n),ifail(k),Z(n,k))

      Call DSYEVX(job,'I',UPLO,n,A,m,VL,VU,1,k,ABSTOL,kk,eval, &
                  Z,n, WORK, LWORK, IWORK, IFAIL, INFO)

      if(kk.ne.k) then
       write(*,*) ' DSYEVX(lapack) provides',kk,' eigenvectors'
       write(*,*) ' when we ordered', k,'  ones'
       Error Stop ' Stop in LAP_DSYEVX'
      end if

      if(INFO.ne.0) write(*,*) 'DSYEVX(lapack) gives INFO = ',INFO

      A(1:n,1:k) = Z(1:n,1:k)

      Deallocate(work,iwork,ifail,Z)

      End Subroutine LAP_DSYEVX


!======================================================================
      Subroutine LAP_DSYGV(job,UPLO,n,m,A,C,eval,info)
!======================================================================
!     Call LAPACK procedure DSYEV to obtain the eigenvalues (eval) and
!     eigenvectors (A) for generalized problem A S = E C S
!     job = 'N' - compute eigenvalues only
!           'V' - compute eigenvalues and eigenvectors
!     UPLO='U'  or 'L' -  used upper or lower triangle matrix
!---------------------------------------------------------------------
      Implicit none
      Character(1), intent(in) :: job,UPLO
      Integer, intent(in) :: n,m
      Integer, intent(out), optional :: info
      Real(8) :: A(m,*),C(m,*)
      Real(8) :: eval(*)

      Integer :: lwork
      Real(8) :: work(3*n)

      lwork = 3*n

      Call DSYGV(1,job,UPLO,n,A,m,C,m,eval,WORK,LWORK,INFO)

      if(INFO.ne.0) write(*,*) ' DSYGV(lapack) gives INFO = ',INFO


      End Subroutine LAP_DSYGV


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
      Implicit none
      Character(1), intent(in) :: job,UPLO
      Integer, intent(in) :: n,m,k
      Integer, intent(out), optional :: info
      Real(8) :: A(m,*),C(m,*)
      Real(8) :: eval(*)

      Real(8), allocatable :: work(:), Z(:,:)
      Integer, allocatable :: iwork(:),ifail(:)

      Integer :: kk, lwork
      Real(8) :: ABSTOL, VL =0.d0, VU = 0.d0
      Real(8), external :: DLAMCH

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

      Deallocate(work,iwork,ifail,Z)

      End Subroutine LAP_DSYGVX


!======================================================================
      Subroutine LAP_INV(n,m,A,info)
!======================================================================
!     Call LAPACK procedures DGETRF and DGETRI to obtain the inverse
!     matrix A^-1 for matrix A
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,m
      Real(8) :: A(m,m)
      Integer, intent(out), optional :: info
      Real(8) :: work(3*n)
      Integer :: ipiv(n), lwork,ierr

      Call DGETRF(n,n,A,m,IPIV,info)

      if(info.ne.0) write(*,*) 'LAP_INV: DGETRF give INFO =',INFO

      lwork = 3*n
      Call DGETRI(n,A,m,IPIV,WORK,LWORK,INFO)

      if(info.ne.0) write(*,*) 'LAP_INV: DGETRI give INFO =',INFO

      End Subroutine LAP_INV


!======================================================================
      Subroutine LAP_INVS(UPLO,n,m,A,info)
!======================================================================
!     Call LAPACK procedures DGETRF and DGETRI to obtain the inverse
!     matrix A^-1 for symmetric matrix
!     UPLO = 'U' - used upper triangle matrix
!            'L' - used lower triangle matrix
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,m
      Character(1), intent(in) :: UPLO
      Real(8) :: A(m,m)
      Integer, intent(out), optional :: info
      Real(8) :: work(3*n)
      Integer :: ipiv(n), lwork,ierr

      lwork = 3*n

      Call DSYTRF(UPLO,n,A,m,IPIV,WORK,LWORK,INFO)

      if(info.ne.0) write(*,*) 'LAP_INV: DGETRF give INFO =',INFO
      if(info.ne.0) Return

      Call DSYTRI(UPLO,n,A,m,IPIV,WORK,INFO)

      if(info.ne.0) write(*,*) 'LAP_INV: DGETRI give INFO =',INFO

      End Subroutine LAP_INVS


!======================================================================
      Subroutine LAP_ZHEEV(job,UPLO,n,m,A,eval,info)
!======================================================================
!     Call LAPACK procedure ZHEEV to obtain the eigenvalues (eval) and
!     eigenvectors (A) for complex Hermitian matrix A(n,n)
!     job = 'V' or 'N' - compute or not the eigenvectors
!     UPLO = 'U' - used upper triangle matrix
!            'L' - used lower triangle matrix
!---------------------------------------------------------------------
      Implicit none
      Character(1), intent(in) :: job,UPLO
      Integer, intent(in) :: n,m
      Integer, intent(out), optional :: info
      COMPLEX(8) :: A(m,m)
      Real(8) :: eval(m)
      Real(8) :: rwork(3*n)
      COMPLEX(8) :: work(2*n)
      Integer :: lwork,lrwork

      lwork = 2*n; lrwork= 3*n

      Call ZHEEV(job,UPLO,n,A,m,eval,WORK,lwork,RWORK,INFO)

      if(INFO.ne.0) write(*,*) ' ZHEEV(lapack) gives INFO = ',INFO

      End Subroutine LAP_ZHEEV
