!======================================================================
      Subroutine Gen_matrix (itype,kpol)
!======================================================================
! ... generate matrix:
!----------------------------------------------------------------------
      Use dbsr_mat
      Use Timer

      Select case(icase)
       Case(0)
          Call TimerStart('Gen_matrix: O_data')
          Call O_data(itype)
          Call TimerStop('Gen_matrix: O_data')
       Case(1)
          Call TimerStart('Gen_matrix: L_data')
          Call L_data(itype)
          Call TimerStop('Gen_matrix: L_data')
       Case(2)
          Call TimerStart('Gen_matrix: R_data 0')
          Call R_data(itype,kpol,0)
          Call TimerStop('Gen_matrix: R_data 0')
          Call TimerStart('Gen_matrix: R_data 1')
          Call R_data(itype,kpol,1)
          Call TimerStop('Gen_matrix: R_data 1')
          Call TimerStart('Gen_matrix: R_data 2')
          Call R_data(itype,kpol,2)
          Call TimerStop('Gen_matrix: R_data 2')
          Call TimerStart('Gen_matrix: R_data 3')
          Call R_data(itype,kpol,3)
          Call TimerStop('Gen_matrix: R_data 3')
       Case(3)
          Call TimerStart('Gen_matrix: S_data')
          Call S_data(itype,kpol)
          Call TimerStop('Gen_matrix: S_data')
      End Select

      End Subroutine Gen_matrix
