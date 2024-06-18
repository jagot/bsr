subroutine timed_section(t0, tend, message)
  implicit none
  Real(8), intent(in) :: t0, tend
  Character(*), intent(in) :: message

  Real(8) :: seconds
  Integer :: hours=0, minutes=0, iseconds

  seconds = tend - t0
  ! write(*,*) "t0:", t0
  ! write(*,*) "tend:", tend
  ! write(*,*) "seconds:", seconds
  if(seconds >= 60) then
     minutes = seconds/60
     seconds = seconds - 60*minutes
  end if
  if(minutes >= 60) then
     hours = minutes/60
     minutes = minutes - 60*hours
  end if
  iseconds = seconds
  seconds = seconds - iseconds

  write(*,'(a,": ",i2,":",i0.2,":",i0.2,f3.2)') message, hours, minutes, iseconds, seconds
end subroutine timed_section

subroutine timed_section_now(t0, message)
  implicit none
  Real(8), intent(in) :: t0
  Character(*), intent(in) :: message

  Real(8) :: tend

  call CPU_time(tend)
  call timed_section(t0, tend, message)
end subroutine timed_section_now
