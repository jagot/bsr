!======================================================================
      Subroutine br_conf_LS1
!======================================================================
!     broadcast data from  modules conf_LS   (without it_state arrays)
!======================================================================
      Use MPI
      Use conf_LS

      Implicit none
      Integer :: myid,ierr

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(mcfg ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ncfg ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kcfg ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lcfg ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) then
       if(allocated(ip_state)) Deallocate(ip_state)
       if(allocated(IC_term )) Deallocate(IC_term )
       if(allocated(ip_orb  )) Deallocate(ip_orb  )
       if(allocated(WC      )) Deallocate(WC      )
       Allocate(ip_state(mcfg),IC_term(mcfg),ip_orb(kcfg),WC(mcfg))
      end if

      Call MPI_BCAST(ip_state,mcfg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(IC_term ,mcfg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ip_orb  ,kcfg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(WC,mcfg,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_conf_LS1

