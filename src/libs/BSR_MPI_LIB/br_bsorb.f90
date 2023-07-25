!======================================================================
      Subroutine br_bsorb
!======================================================================
!     broadcast data from  module spline_orbitals
!======================================================================
      Use MPI
      Use spline_orbitals
      Use spline_param, only: ns

      Implicit none
      Integer :: myid,ierr
      
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(nbf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(nbf.eq.0) Return

      if(myid.ne.0) then
       if(Allocated(nbs)) & 
        Deallocate (nbs,lbs,kbs,mbs,iech,ebs,PBS,QBS)
       mbf = nbf
       Allocate(nbs(mbf),lbs(mbf),kbs(mbf),ebs(mbf),mbs(1:mbf), &
                iech(1:mbf),PBS(1:ns,1:mbf),QBS(1:ns,1:mbf))
      end if

      Call MPI_BCAST(nbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(iech,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ebs,nbf*4,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(PBS,nbf*ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(QBS,nbf*ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


      Call MPI_BCAST(nv_ch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(nv_ch.eq.0) Return

      if(myid.ne.0) then
       if(Allocated(i_ch)) Deallocate (i_ch)
       if(Allocated(j_ch)) Deallocate (j_ch)
       if(Allocated(V_ch)) Deallocate (V_ch)
       Allocate(i_ch(nv_ch),j_ch(nv_ch),V_ch(ns,nv_ch))
      end if

      Call MPI_BCAST(i_ch,nv_ch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(j_ch,nv_ch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(V_ch,nv_ch*ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_bsorb

