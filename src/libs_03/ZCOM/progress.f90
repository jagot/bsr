module ProgressMeter
  type Progress
     integer :: N ! Number of steps
     real(8) :: tstart
     real(8) :: dt_print
     real(8) :: tlast_print
  end type Progress
contains
  subroutine init_progress(p, N)
    implicit none
    type(Progress), intent(inout) :: p
    integer, intent(in) :: N

    p%N = N
    p%dt_print = 0.1
    call CPU_time(p%tstart)
    p%tlast_print = p%tstart
    write(*,*)
  end subroutine init_progress

  subroutine carriage_return
    write(*,'(a)', advance='no') achar(13)
  end subroutine carriage_return

  subroutine clear_line(n)
    implicit none
    integer, intent(in) :: n
    integer :: i

    call carriage_return
    do i=1,n
       write(*,'(a)', advance='no') " "
    end do
    call carriage_return
  end subroutine clear_line

  subroutine print_progress(p, i)
    use Timer
    implicit none
    type(Progress), intent(inout) :: p
    integer, intent(in) :: i
    real(8) :: tnow, elapsed, eta

    call CPU_time(tnow)
    elapsed = tnow - p%tstart
    if(tnow - p%tlast_print < p%dt_print) return

    eta = (elapsed*p%N)/i

    call clear_line(100)
    write(*,'("Progress: ", f5.2, " %")', advance='no') (100.0*i)/p%N
    write(*,'(a)', advance='no') " "
    call print_time(elapsed, "Elapsed")
    write(*,'(a)', advance='no') " "
    call print_time(eta, "ETA")

    p%tlast_print = tnow
  end subroutine print_progress
end module ProgressMeter
