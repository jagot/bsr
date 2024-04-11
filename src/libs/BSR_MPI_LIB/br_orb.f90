!======================================================================
      Subroutine br_orb
!======================================================================
!     broadcast data from  module spline_orbitals
!======================================================================
      Use MPI
      Use orb_LS

      Implicit none
      Integer :: myid,ierr
    
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(nwf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      mwf = nwf
      if(nwf.eq.0) Return

      if(myid.ne.0) then
       if(Allocated(NEF)) Deallocate (NEF,LEF,KEF,IEF,ELF,IORT)
       Allocate(NEF(mwf),LEF(mwf),KEF(mwf),IEF(mwf),ELF(mwf),IORT(mwf,mwf))
      end if

      Call MPI_BCAST(NEF ,mwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(LEF ,mwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(KEF ,mwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IEF ,mwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ELF,mwf*4,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(IORT,mwf*mwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_orb

