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

module SimpleProgress
  use accuracy
  use Timer, only: get_real_time, print_time

  implicit none
  private

  public simple_progress, simple_progress_header, simple_progress_start

  interface simple_progress
     procedure simple_progress_full, simple_progress_default_advance, simple_progress_advance_only, simple_progress_min
  end interface simple_progress

  interface simple_progress_header
     procedure simple_progress_header_full, simple_progress_header_advance
  end interface simple_progress_header

  interface simple_progress_start
     procedure simple_progress_start_full, simple_progress_start_advance
  end interface simple_progress_start

contains
  subroutine seconds2hhmmss(hours, minutes, seconds)
    integer(ik), intent(out) :: hours, minutes
    real(rk), intent(inout) :: seconds
    hours = 0
    minutes = 0
    if(seconds > 60) then
       minutes = int(seconds/60)
       seconds = seconds - 60*minutes
    end if
    if(minutes > 60) then
       hours = minutes/60
       minutes = minutes - 60*hours
    end if
  end subroutine seconds2hhmmss

  subroutine simple_progress_full(i, num_steps, starting_time, perf_multiplier, advance)
    integer(ik), intent(in) :: i, num_steps
    real(rk), intent(in) :: starting_time
    ! perf_multiplier can be e.g. the number of grid points, then
    ! the quantity that is printed corresponds to space--time grid
    ! points per second; this is inspired by Ken Schafer.
    real(rk), intent(in) :: perf_multiplier
    logical, intent(in) :: advance

    real(rk) perf_multiplier_, elapsed, seconds_per_step, eta
    integer(ik) el_hours, el_minutes, eta_hours, eta_minutes
    Character(3) :: advance_ = 'yes'

    if(.not.advance) then
       advance_ = 'no'
    end if

    elapsed = get_real_time() - starting_time

    seconds_per_step = elapsed / i
    eta = max(num_steps*seconds_per_step - elapsed, 0.0_rk)

    call seconds2hhmmss(el_hours, el_minutes, elapsed)
    call seconds2hhmmss(eta_hours, eta_minutes, eta)

    write(out, '(i10, 1i4.2,":",1i2.2,":",1i2.2, 1g26.16, " Hz", 1i4.2,":",1i2.2,":",1i2.2)', &
         advance=trim(advance_)) &
         i, el_hours, el_minutes, int(elapsed), &
         perf_multiplier/seconds_per_step, &
         eta_hours, eta_minutes, int(eta)
  end subroutine simple_progress_full

  subroutine simple_progress_default_advance(i, num_steps, starting_time, perf_multiplier)
    integer(ik), intent(in) :: i, num_steps
    real(rk), intent(in) :: starting_time
    real(rk), intent(in) :: perf_multiplier
    call simple_progress_full(i, num_steps, starting_time, perf_multiplier, .true.)
  end subroutine simple_progress_default_advance

  subroutine simple_progress_advance_only(i, num_steps, starting_time, advance)
    integer(ik), intent(in) :: i, num_steps
    real(rk), intent(in) :: starting_time
    logical, intent(in) :: advance
    call simple_progress_full(i, num_steps, starting_time, 1.0_rk, advance)
  end subroutine simple_progress_advance_only

  subroutine simple_progress_min(i, num_steps, starting_time)
    integer(ik), intent(in) :: i, num_steps
    real(rk), intent(in) :: starting_time
    call simple_progress_default_advance(i, num_steps, starting_time, 1.0_rk)
  end subroutine simple_progress_min

  subroutine simple_progress_header_full(advance)
    logical, intent(in) :: advance
    Character(3) :: advance_ = 'yes'
    if(.not.advance) then
       advance_ = 'no'
    end if
    write(out, '(a10, a10, a29, a10)', advance=trim(advance_)) &
         "Step", "Elapsed", "Performance", "ETA"
  end subroutine simple_progress_header_full

  subroutine simple_progress_header_advance()
    call simple_progress_header_full(.true.)
  end subroutine simple_progress_header_advance

  function simple_progress_start_full(advance) result(now)
    logical, intent(in) :: advance
    real(rk) :: now

    call simple_progress_header(advance)
    now = get_real_time()
  end function simple_progress_start_full

  function simple_progress_start_advance() result(now)
    real(rk) :: now
    now = simple_progress_start_full(.true.)
  end function simple_progress_start_advance
end module SimpleProgress
