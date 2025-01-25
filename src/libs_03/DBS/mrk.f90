subroutine mrk_common_gen_array(ksi, ksj, rkd1, rkd2, rkd3, rkd4)
  Use DBS_grid
  Use DBS_moments, only: rkd
  Use DBS_integrals
  Use Timer

  Implicit none
  Integer, intent(in) :: ksi, ksj
  Real(8), intent(in) :: rkd1(1:ks*ks,1:nv), rkd2(1:ks*ks,1:nv), rkd3(1:ks*ks,1:nv), rkd4(1:ks*ks,1:nv)

  Integer :: lii_max, lij_max
  Integer :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
  Real(8) :: c

  call TimerStart('mrk: generate array')
  rkb_tmp = 0.d0
  ! Compute maximum linear indices
  lij_max = ksj*(ksj+1)/2
  lii_max = ksi*(ksi+1)/2
  do jj = 1,lij_max
     do ii = 1,lii_max
        !$omp parallel do default(none) &
        !$omp private(jv,iv, jhp,ihp, jh,ih, jp,ip, j,i, c) &
        !$omp shared(jj,ii, ksj,ksi, nv, rkd1,rkd2,rkd3,rkd4,rkd, rkb_tmp)
        do jv = 1,nv
           call mrk_compute_from_linear_index(ksj, jj, jhp, jh, jp)
           call mrk_compute_from_linear_index(ksi, ii, ihp, ih, ip)
           j = jv + jh - 1

           ! Sub-diagonal
           do iv=1,jv-1
              i = iv + ih - 1
              c = rkd1(ii,iv)*rkd2(jj,jv)
              rkb_tmp(i,j,ip,jp) = rkb_tmp(i,j,ip,jp) +  c
           end do

           ! Super-diagonal
           do iv=jv+1,nv
              i = iv + ih - 1
              c = rkd3(jj,jv)*rkd4(ii,iv)
              rkb_tmp(i,j,ip,jp) = rkb_tmp(i,j,ip,jp) +  c
           end do

           ! Diagonal
           do iv = jv,jv
              i = iv + ih - 1
              rkb_tmp(i,j,ip,jp) = rkb_tmp(i,j,ip,jp) + rkd(ii,jj,iv)
           end do
        end do
        !$omp end parallel do
     end do
  end do
  rkb(1:ns,1:ns,1:ks,1:ks) = rkb_tmp(1:ns,1:ns,1:ks,1:ks)
  call TimerStop('mrk: generate array')

contains
  subroutine mrk_compute_indices(ks, iv, ih, ihp, i, ip, ii)
    Implicit none
    Integer, intent(in) :: ks, iv, ih, ihp
    Integer, intent(out) :: i, ip, ii

    i  = iv  + ih - 1
    ip = ihp - ih + 1

    ii = ihp - (ih*(-1 + ih - 2*ks))/2 - ks
  end subroutine mrk_compute_indices

  subroutine mrk_compute_indices_small(ks, ih, ihp, ip, ii)
    Implicit none
    Integer, intent(in) :: ks, ih, ihp
    Integer, intent(out) :: ip, ii

    ip = ihp - ih + 1
    ii = ihp - (ih*(-1 + ih - 2*ks))/2 - ks
  end subroutine mrk_compute_indices_small

  subroutine mrk_compute_from_linear_index(n, li, i, j, ip)
    Implicit none
    Integer, intent(in) :: n, li
    Integer, intent(out) :: i, j, ip
    call linear_to_triangle_indices(n, li, i, j)
    ip = i - j + 1
  end subroutine mrk_compute_from_linear_index

  function triangle_to_linear_index(k, i, j) result(li)
    Implicit none
    Integer, intent(in) :: k, i, j
    Integer :: li

    li = i - (j*(-1 + j - 2*k))/2 - k
  end function triangle_to_linear_index

! https://atrebas.github.io/post/2021-01-17-index_to_lower_triangular_subscripts/
! Cf. Functions for col-wise numbering, diagonal included
  subroutine linear_to_triangle_indices(n, k, i, j)
    Use accuracy, only: rk
    Implicit none
    Integer, intent(in) :: n ! Size of the lower-triangular matrix
    Integer, intent(in) :: k ! Linear index
    Integer, intent(out) :: i, j ! Row and column indices
    Integer :: kp, p

    kp = n*(n+1)/2 - k
    p = floor((sqrt(1 + 8._rk*kp)-1)/2)
    i = n - (kp - p*(p+1)/2)
    j = n - p

    ! write(*, '(6i5)') n, k, kp, p, i, j

  end subroutine linear_to_triangle_indices

  subroutine test_triangular_indices(ksj)
    Implicit none
    Integer, intent(in) :: ksj
    Integer :: jj, jh, jhp, i, j
    write(*,'("ksj = ",i0)') ksj
    jj = 0
    do jh = 1,ksj
       do jhp=jh,ksj
          jj = jj + 1
          call linear_to_triangle_indices(ksj, jj, i, j)
          write(*,'("jh = ",i0,", jhp = ",i0,", jj =",2(xi3),", i = ",i0,", j = ",i0)') &
               jh, jhp, jj, &
               triangle_to_linear_index(ksj, jhp, jh), &
               i, j
          if(j /= jh) error stop "Linear to triangular index failed"
          if(i /= jhp) error stop "Linear to triangular index failed"
       end do
       write(*,*)
    end do
    error stop "Booh"
  end subroutine test_triangular_indices
end subroutine mrk_common_gen_array
