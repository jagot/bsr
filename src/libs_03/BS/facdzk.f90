!======================================================================
      Subroutine FACDZK (ktx,k,dzk,ipvtz)
!======================================================================
!
!     Sets up and factorizes the matrix
!
!         D(i,j) + k/r
!
!     for solving the following differential equation by the
!     spline-Galerkin method:
!
!         (d + k/r) zk(r) = f(r)
!
!     Calls:  dgbtrf (LAPACK) or  dgbfa (LINPACK)
!
!     on entry
!     --------
!     k      define equation
!
!     on exit
!     -------
!     dzk  factorized array for the differential with operator
!          (d+k/r).  dzk is banded with width 3*ks-2.
!          The first ks-1 rows are zeros, and the
!          following k rows are the elements above the diagnal,
!          the last ks-1 rows are the elements below the diagonal of
!          the original arrays.
!     iptvz  an integer vector of pivot indices.
!--------------------------------------------------------------------

    USE spline_param
    USE spline_galerkin

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ktx,k
    REAL(KIND=8), DIMENSION(ktx,ns), INTENT(OUT) :: dzk
    INTEGER, DIMENSION(ns), INTENT(OUT) :: ipvtz

    ! ..  local variables

    INTEGER :: i, j, ier


! .. clear the first empty rows of dzk array  ...

     do i=1,ks-1
      do j=1,ns
        dzk(i,j)=0.d0
      end do
     end do

! .. lower portion ...  (last rows of dyk)

     do j=1,ks
       do i=ks-j+1,ns
         dzk(3*ks-1-j,i-ks+j) = db1(i,j) + k*rm1(i,j)
       end do
     end do

! .. upper portion ...

     do j=2,ks
       do i=1,ns-j+1
         dzk(2*ks-j,i+j-1) = - db1(i+j-1,ks-j+1) + k*rm1(i+j-1,ks-j+1)
       end do
     end do

! .. preparation for boundary conditions

      do i=1,ks
       j=2*ks-i
       dzk(j,i)=0.d0                  ! first equation
      end do                          ! for zero value in origin
      dzk(2*ks-1,1)=1.d0

      if(k.eq.-1) then
       do i=1,ks
        dzk(2*ks+1-i,i)=0.d0          ! second equation (for Nk integrals)
       end do                         ! for zero derivative in origin
       dzk(2*ks-1,2)=1.d0
      end if

      if(k.le.-3) then
       do i=1,ks
        dzk(3*ks-1-i,ns-ks+i)=0.d0    ! last equation (for Vk integrals)
       end do                         ! for boudary condition at the end
       dzk(2*ks-1,ns) = 1.d0
      end if

! .. factorize dzk ...

!      CALL dgbfa(dzk,ktx,ns,ks-1,ks-1,ipvtz,ier)
!      IF (ier /= 0) STOP 'FACDZK: dgbfa (LINPACK) failed'

      Call DGBTRF (ns,ns,ks-1,ks-1,dzk,3*ks-2,ipvtz,ier)
      if (ier .ne. 0) stop 'FACDZK: dgbtrf (LAPACK) failed'

    End Subroutine FACDZK
