module Timer
contains
  subroutine print_time(seconds, message)
    use, intrinsic :: iso_fortran_env, only : stdout=>output_unit

    implicit none
    Real(8), intent(in) :: seconds
    Character(*), intent(in) :: message

    Real(8) :: sec
    Integer :: hours, minutes, iseconds

    hours = 0
    minutes = 0
    iseconds = 0

    sec = seconds
    if(sec >= 60) then
       minutes = sec/60
       sec = sec - 60*minutes
    end if
    if(minutes >= 60) then
       hours = minutes/60
       minutes = minutes - 60*hours
    end if
    iseconds = sec
    sec = sec - iseconds

    write(stdout,'(a,": ",i2,":",i0.2,":",i0.2,f3.2)', advance='no') message, hours, minutes, iseconds, sec
  end subroutine print_time

  subroutine timed_section(t0, tend, message)
    use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
    implicit none
    Real(8), intent(in) :: t0, tend
    Character(*), intent(in) :: message
    call print_time(tend - t0, message)
    write(stdout,*)
    call flush(stdout)
  end subroutine timed_section

  subroutine timed_section_now(t0, message)
    implicit none
    Real(8), intent(in) :: t0
    Character(*), intent(in) :: message

    Real(8) :: tend

    call CPU_time(tend)
    call timed_section(t0, tend, message)
  end subroutine timed_section_now
end module Timer
