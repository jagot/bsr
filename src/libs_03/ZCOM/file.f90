! Taken from https://fortran-lang.discourse.group/t/do-opened-files-close-automatically-when-they-go-out-of-scope/4065/4
module file_helper
  type :: file
     integer :: unit = -1
   contains
     final :: close_file
  end type file
contains
  subroutine close_file(this)
    type(file), intent(inout) :: this

    if (this%unit /= -1) then
       close(this%unit)
    end if

  end subroutine close_file
end module file_helper
