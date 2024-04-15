!======================================================================
      Subroutine Read_cfgu(nu)
!======================================================================
!     Read the configuration list and auxiliary arrays from unformatted
!     c-file
!----------------------------------------------------------------------
      Use symc_list_LS
      Use symt_list_LS
      Use conf_LS
      Use orb_LS

      Implicit none
      Integer, intent(in) :: nu
      Integer :: i

      Call Read_orb_LS(nu,nclosd)
      Call Read_symc_LS(nu)
      Call Read_symt_LS(nu)
      Call Read_conf_LS(nu)

      if(allocated(IT_sort )) Deallocate(IT_sort)
      if(allocated(IC_term1)) Deallocate(IC_term1)
      if(allocated(IC_term2)) Deallocate(IC_term2)
      Allocate(IT_sort(nsymt),IC_term1(nsymc),IC_term2(nsymc))

      read(nu) (IT_sort(i),i=1,nsymt)
      read(nu) (IC_term1(i),i=1,nsymc)
      read(nu) (IC_term2(i),i=1,nsymc)

      if(allocated(IP_stat  )) Deallocate(IP_stat  )
      if(allocated(IT_state1)) Deallocate(IT_state1)
      if(allocated(IT_state2)) Deallocate(IT_state2)
      Allocate(IP_stat(ncfg),IT_state1(nsymt),IT_state2(nsymt))

      read(nu) (IP_stat(i),i=1,ncfg)
      read(nu) (IT_state1(i),i=1,nsymt)
      read(nu) (IT_state2(i),i=1,nsymt)

      End Subroutine Read_cfgu



