module BLAS
  implicit none
  private

  public general_rank_update, &
       blas_set_num_threads, blas_get_num_threads

  interface general_rank_update
     module procedure sgeneral_rank_update
     module procedure dgeneral_rank_update
     module procedure cgeneral_rank_update
     module procedure zgeneral_rank_update
  end interface general_rank_update

  interface
#ifndef BLA_VENDOR
#error "You must set BLA_VENDOR as a compile definition"
#elif BLA_VENDOR == OpenBLAS
     subroutine blas_set_num_threads(n) bind(C, name='openblas_set_num_threads')
       use ISO_C_BINDING
       integer(C_INT), intent(in), value :: n
     end subroutine blas_set_num_threads

     integer(C_INT) function blas_get_num_threads() bind(C, name='openblas_get_num_threads')
       use ISO_C_BINDING
     end function blas_get_num_threads
#else
#error "Implementation of blas_{s,g}et_num_threads missing for chosen BLA_VENDOR"
#endif
  end interface
contains
  ! * BLAS 2
  ! ** General rank update
  subroutine sgeneral_rank_update(alpha, x, y, A)
    real(4), intent(in) :: alpha, x(:), y(:)
    real(4), intent(inout) :: A(:,:)
    integer :: m, n
    m = size(A, dim=1)
    n = size(A, dim=2)
    call sger(m, n, alpha, x, 1, y, 1, A, m)
  end subroutine sgeneral_rank_update

  subroutine dgeneral_rank_update(alpha, x, y, A)
    real(8), intent(in) :: alpha, x(:), y(:)
    real(8), intent(inout) :: A(:,:)
    integer :: m, n
    m = size(A, dim=1)
    n = size(A, dim=2)
    call dger(m, n, alpha, x, 1, y, 1, A, m)
  end subroutine dgeneral_rank_update

  subroutine cgeneral_rank_update(alpha, x, y, A)
    complex(8), intent(in) :: alpha, x(:), y(:)
    complex(8), intent(inout) :: A(:,:)
    integer :: m, n
    m = size(A, dim=1)
    n = size(A, dim=2)
    call cgerc(m, n, alpha, x, 1, y, 1, A, m)
  end subroutine cgeneral_rank_update

  subroutine zgeneral_rank_update(alpha, x, y, A)
    complex(16), intent(in) :: alpha, x(:), y(:)
    complex(16), intent(inout) :: A(:,:)
    integer :: m, n
    m = size(A, dim=1)
    n = size(A, dim=2)
    call zgerc(m, n, alpha, x, 1, y, 1, A, m)
  end subroutine zgeneral_rank_update
end module BLAS
