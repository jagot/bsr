module LAPACK
  implicit none
  private

  public herm_eig_mrrr
  interface herm_eig_mrrr
     module procedure dherm_eig_mrrr
     module procedure dherm_eig_mrrr_all
     module procedure dherm_eig_mrrr_int_range
     module procedure dherm_eig_mrrr_real_interval
  end interface herm_eig_mrrr

contains
  function dherm_eig_mrrr(job,range,UPLO,A,eval,IL,IU,VL,VU) result(M)
    Implicit none
    Character(1), intent(in) :: job,range,UPLO
    Real(8), intent(inout) :: A(:,:)
    Real(8), intent(inout) :: eval(:)
    Integer, intent(in) :: IL,IU
    Real(8), intent(in) :: VL, VU
    Integer :: M

    Integer :: n

    Real(8), allocatable :: work(:), Z(:,:)
    Integer, allocatable :: iwork(:),isuppz(:)
    Integer :: lwork, liwork, info
    Real(8) :: ABSTOL
    Real(8), external :: DLAMCH
    Character(256) :: error_msg

    n = size(A,1)

    if(range == 'I'.and.(IL<1.or.IU<IL)) then
       write(error_msg, '("IL,IU = ",i0,","i0,"; 1 <= IL <= IU required")') IL,IU
       error stop trim(error_msg)
    end if
    if(range == 'V'.and.(VL>=VU)) then
       write(error_msg, '("VL,VU = ",e26.16,","e26.16,"; VL < VU required")') VL,VU
       error stop trim(error_msg)
    end if

    ABSTOL = 4*DLAMCH('S')
    lwork = 26*n
    liwork = 10*n
    M = n
    if(range == 'I') then
       M = IU-IL+1
    end if
    Allocate(work(lwork),iwork(liwork),Z(n,M),isuppz(2*M))

    Call DSYEVR(job,range,UPLO,n,A,n,VL,VU,IL,IU,ABSTOL,M,eval, &
         Z,n, isuppz, WORK, LWORK, IWORK, liwork, INFO)
    if(info /= 0) then
       write(error_msg, '("DSYEVR error code INFO = ",i0)') info
       error stop trim(error_msg)
    end if

    if(range == 'I') then
       if(M.ne.(IU-IL+1)) then
          write(*,*) ' DSYEVR(lapack) provides',M,' eigenvectors'
          write(*,*) ' when we ordered', (IU-IL+1),'  ones'
          Error Stop ' Stop in LAP_DSYEVR'
       end if
    end if

    A(1:n,1:M) = Z(1:n,1:M)

    Deallocate(work,iwork,Z,isuppz)
  end function dherm_eig_mrrr

  function dherm_eig_mrrr_all(job,UPLO,A,eval) result(M)
    Implicit none
    Character(1), intent(in) :: job,UPLO
    Real(8), intent(inout) :: A(:,:)
    Real(8), intent(inout) :: eval(:)
    Integer :: M

    M = dherm_eig_mrrr(job,'A',UPLO,A,eval,0,0,0.d0,0.d0)
  end function dherm_eig_mrrr_all

  function dherm_eig_mrrr_int_range(job,UPLO,A,eval,IL,IU) result(M)
    Implicit none
    Character(1), intent(in) :: job,UPLO
    Real(8), intent(inout) :: A(:,:)
    Real(8), intent(inout) :: eval(:)
    Integer, intent(in) :: IL,IU
    Integer :: M

    M = dherm_eig_mrrr(job,'I',UPLO,A,eval,IL,IU,0.d0,0.d0)
  end function dherm_eig_mrrr_int_range

  function dherm_eig_mrrr_real_interval(job,UPLO,A,eval,VL,VU) result(M)
    Implicit none
    Character(1), intent(in) :: job,UPLO
    Real(8), intent(inout) :: A(:,:)
    Real(8), intent(inout) :: eval(:)
    Real(8), intent(in) :: VL, VU
    Integer :: M

    M = dherm_eig_mrrr(job,'V',UPLO,A,eval,0,0,VL,VU)
  end function dherm_eig_mrrr_real_interval
end module LAPACK
