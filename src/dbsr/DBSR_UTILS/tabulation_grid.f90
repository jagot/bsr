module grid_tools
  implicit none
contains
  subroutine tabulation_grid(r, rmax, distribution)
    implicit none

    real(8), intent(out) :: r(:)
    real(8), intent(in) :: rmax
    integer, intent(in) :: distribution

    real(8) :: dr, tr
    integer :: nr, i

    nr = size(r, 1)

    dr = 1.0/(nr-1)

    do i=1,nr
       tr = (i-1)*dr
       if(distribution==1) then
       elseif(distribution==2) then
          tr = tr*tr
       elseif(distribution==3) then
          tr = (exp(tr)-1)/(exp(1.0)-1)
       else
          error stop 'Unknown distribution'
       end if
       r(i) = tr * rmax
    end do
  end subroutine tabulation_grid

  subroutine print_distribution_help
    implicit none

    write(*,*) "  distribution is an integer deciding the spacing of the radial grid points:"
    write(*,*) ""
    write(*,*) "  r = f(t)*rmax"
    write(*,*) ""
    write(*,*) "    1: f(t) = t"
    write(*,*) "    2: f(t) = t^2"
    write(*,*) "    3: f(t) = (e^t-1)/(e-1)"
    write(*,*) ""
    write(*,*) "  where t varies from 0 to 1"
  end subroutine print_distribution_help
end module grid_tools
