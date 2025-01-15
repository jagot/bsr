subroutine assert_different_files(src, dst, action, preposition)
  Implicit none

  Character(len=*) :: src, dst, action, preposition
  Character(256) error_msg

  ! This does not cover the possibility of one file being a symlink to
  ! the other.
  if(trim(src) == trim(dst)) then
     write(error_msg, '("Cannot ",a," ",a," ",a," identical file ",a)') &
          trim(action), trim(src), trim(preposition), trim(dst)
     error stop error_msg
  end if

end subroutine assert_different_files

subroutine copy_file_data(us, ud)
  Implicit none

  Integer :: us, ud

  ! Integer(1) is actually platform-dependent, not necessarily a
  ! byte. For most compilers it does, but no guarantees.
  Integer(1) :: buffer(16*1024)
  Integer :: count, stat, pos1, pos2

  do
     ! inquire returns pos in "storage units", not necessarily
     ! bytes. However, for access='stream', a storage unit is
     ! effectively a byte.
     inquire(us,pos=pos1)
     read(us,iostat=stat) buffer
     inquire(us,pos=pos2)
     count = pos2 - pos1
     ! write(*,*) count, stat
     write(ud) buffer(1:count)
     if(stat /= 0) exit
  end do

end subroutine copy_file_data

! This subroutine merges the contents of two files dst and src and
! stores the results in the file dst, such that its contents appears
! before those of src. This is more portable than using System calls
! to cat, which only works under Unix.
subroutine merge_files(src, dst)
  Implicit none

  Character(len=*) :: src, dst

  Integer :: us, ud

  call assert_different_files(src, dst, 'merge', 'with')

  ! write(*,'(a,": ",a," >> ",a)') 'We are trying to merge files here', &
  !      trim(src), trim(dst)

  open(newunit=ud, file=trim(dst), access='stream', status='old', &
       form='unformatted', action='readwrite', position='append')
  open(newunit=us, file=trim(src), access='stream', status='old', &
       form='unformatted', action='read', position='rewind')

  call copy_file_data(us, ud)

  close(ud)
  close(us)
end subroutine merge_files

subroutine move_file(src, dst)
  Implicit none

  Character(len=*) :: dst, src

  Integer :: us, ud

  call assert_different_files(src, dst, 'move', 'to')

  ! write(*,'(a,": ",a," -> ",a)') 'We are trying to move a file here', &
  !      trim(src), trim(dst)

  open(newunit=ud, file=trim(dst), access='stream', status='replace', &
       form='unformatted', action='write', position='rewind')
  open(newunit=us, file=trim(src), access='stream', status='old', &
       form='unformatted', action='read', position='rewind')

  call copy_file_data(us, ud)

  close(ud)
  close(us)
end subroutine move_file

subroutine delete_file(file)
  Implicit none
  Character(len=*) file
  Integer :: uf

  ! write(*,'(a,": ",a)') 'We are trying to delete a file here', &
  !      trim(file)

  open(newunit=uf, file=trim(file), status='old')
  close(uf,status='delete')
end subroutine delete_file
