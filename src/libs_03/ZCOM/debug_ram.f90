!
!  Debugging routine: Report total VM size. Linux-only!
!
module debug_ram
  use accuracy
  public trace_ram
  public get_ram_size
  public rcsid_debug_ram
  !
  !
  character(len=clen), save :: rcsid_debug_ram = "$Id: debug_ram.f90,v 1.6 2022/11/17 06:26:05 ps Exp $"
  real(xrk), save           :: old_vmsize         = -1._rk
  !
  contains

  subroutine trace_ram(tag)
    character(len=*), intent(in) :: tag
    real(xrk)                    :: vmsize
    !
    vmsize = get_ram_size()
    if (vmsize/=old_vmsize) then
      write (out,"('#DBG: VM size changed from ',f0.0,' to ',f0.0,' tag= ',a)") &
             old_vmsize, vmsize, trim(tag)
      call flush_wrapper(out)
      old_vmsize = vmsize
    end if
  end subroutine trace_ram
  !
  function get_ram_size() result(vmsize)
    real(xrk)           :: vmsize  ! Total VM size we currently use, in kilobytes
    !
    integer(ik)         :: iu, ios
    integer(hik)        :: part_size
    character(len=clen) :: buf
    !
    iu = -1
    vmsize = 0
    error_block: do
      open (newunit=iu,file="/proc/self/smaps",form="formatted",action="read",iostat=ios)
      if (ios/=0) exit error_block
      read_loop: do
        read (iu,"(a)",iostat=ios) buf
        if (ios/=0) exit read_loop
        if (buf(1:12)/="Referenced: ") cycle read_loop
        read (buf(13:),*) part_size
        vmsize = vmsize + real(part_size,kind=xrk)
      end do read_loop
      close(iu)
      return
    end do error_block
    ! write (out,"('get_ram_size: Error ',i0)") ios
    ! call flush_wrapper(out)
    if (iu/=-1) close (iu)
  end function get_ram_size
end module debug_ram
