!=====================================================================
      Subroutine Sort_states
!---------------------------------------------------------------------
! ... sorting states according their term symmetries
!     (define IT_state1 and IT_state2 arrays in the conf_LS module):
!---------------------------------------------------------------------
      Use symt_list_LS
      Use conf_LS

      Implicit none
      Integer :: i,it

      if(allocated(IP_stat  )) Deallocate(IP_stat  )
      if(allocated(IT_state1)) Deallocate(IT_state1)
      if(allocated(IT_state2)) Deallocate(IT_state2)
      Allocate(IP_stat(ncfg),IT_state1(nsymt),IT_state2(nsymt))
      IT_state1=0; IT_state2=0

      Call SORT_IT (ncfg,IC_term,IP_stat)

      Do i=1,ncfg
       it=IC_term(IP_stat(i))
       if(IT_state1(it).eq.0) IT_state1(it)=i; IT_state2(it)=i 
      End do

      End Subroutine Sort_states

