!=====================================================================
      Subroutine SUB1_HD
!=====================================================================
!     calculations for given partial wave
!---------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Integer :: num_sol ! Number of solutions found

! ... read channel information:

      Call Read_channel_jj(nut,klsp)

! ... check and print main parameters and file:

      Call Check_dbsr_mat

! ... diagonalize the matrix and get inner-region solutions:

      Call Diag_mat(num_sol)

! ... output of solutions and find the surface amplitudes:

      if(itype.ge. 0)  Call Rsol_out(num_sol)

! ... output of standard H.nnn file:

      if(itype.ge. 0 .and. iiexp.eq.0)  Call H_out(num_sol)
      if(itype.ge. 0 .and. iiexp.gt.0)  Call H_out_exp(num_sol)

! ... output of bound states in bound.nnn:

      if(itype.eq.-1)  Call B_out(num_sol)

      End Subroutine SUB1_HD
