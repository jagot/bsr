!======================================================================
      Subroutine br_ovl
!======================================================================
!     broadcast data from  module spline_orbitals
!======================================================================
      Use MPI
      Use orb_overlaps

      Implicit none
      Integer :: myid,ierr
      
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(norb, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(max_l,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nobs ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      mem_orb_overlaps = 2*nobs + 4*norb + 3*max_l

      if(norb.eq.0) Return

      if(myid.ne.0) then
       if(allocated(ipl)) Deallocate(ipl,jpl,lorb,chan,ip_l,jp_l,ip_ovl,Cobs)
       Allocate(lorb(norb),chan(norb),ipl(norb),jpl(norb))
       Allocate(ip_l(0:max_l),jp_l(0:max_l),ip_ovl(0:max_l))
       Allocate(Cobs(nobs))
      end if

      Call MPI_BCAST(lorb,norb,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(chan,norb,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ipl ,norb,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jpl ,norb,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ip_l  (0:max_l),max_l+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jp_l  (0:max_l),max_l+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ip_ovl(0:max_l),max_l+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(Cobs,nobs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_ovl

