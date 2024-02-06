!======================================================================
      Subroutine send_res(ic,jc)
!======================================================================
      Use MPI

      Use mult_par,  only: myid, ierr, JT_oper
      Use term_exp,  only: kt1,kt2, IP_kt1,IP_kt2
      Use ndet_list, only: ndet,ldet,KPD,IPD,NPD
      Use ndef_list, only: ndef,ldef,KPF,IPF,NPF
      Use coef_list, only: ntrm,ncoef,idfc,intc,coef

      Implicit none
      Integer :: ic,jc

      Call MPI_SEND(ncoef,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(ic,   1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(jc,   1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(kt1,  1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(kt2,  1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(ntrm, 1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(JT_oper,ntrm,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

! ... term indexes:

      Call MPI_SEND(IP_kt1,kt1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IP_kt2,kt2,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      if(ncoef.eq.0) Return

! ... coef.s:

      Call MPI_SEND(idfc,ncoef,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(intc,ncoef,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(coef,ntrm*ncoef,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)

! ... dets:

      Call MPI_SEND(ndet,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      if(ndet.eq.0) Return
      Call MPI_SEND(ldet,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(KPD,ndet,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IPD,ndet,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(NPD,ldet,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

! ... defs:

      Call MPI_SEND(ndef,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      if(ndef.eq.0) Return
      Call MPI_SEND(ldef,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      Call MPI_SEND(KPf,ndef,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(IPf,ndef,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(NPf,ldef,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      End Subroutine send_res
