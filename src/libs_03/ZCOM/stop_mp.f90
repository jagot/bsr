!======================================================================
    Subroutine Stop_mpi (pri, info, msg)
!======================================================================
!    Use MPI

    Integer, intent(in)      :: pri, info
    Character(*), intent(in) :: msg

!    Integer :: myid, ierr
!    Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
!    write (pri,*) 'myid = ', myid
!    Call MPI_ABORT(MPI_COMM_WORLD,1,ierr)

    write (pri,'(a)')  trim(msg)
    if(info.ne.0) write (pri,'(/a,i8)') 'info = ',info
    Stop

    End Subroutine Stop_mpi
