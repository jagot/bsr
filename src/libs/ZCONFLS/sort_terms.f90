!======================================================================
      Subroutine Sort_terms
!----------------------------------------------------------------------
! ... sorting the terms according to configuration symmetries,
! ... (define IC_term1 and IC_term2 arrays in symc_list_LS module)
!----------------------------------------------------------------------
      Use  symc_list_LS
      Use  symt_list_LS

      Implicit none
      Integer :: it, ic 

      if(allocated(IT_sort )) Deallocate(IT_sort)
                              Allocate(IT_sort(nsymt))
      if(allocated(IC_term1)) Deallocate(IC_term1)
                              Allocate(IC_term1(nsymc))
      if(allocated(IC_term2)) Deallocate(IC_term2)
                              Allocate(IC_term2(nsymc))

      Call Sort_IT(nsymt,IT_conf,IT_sort)

      IC_term1=0
      Do it=1,nsymt
       ic=IT_conf(IT_sort(it))
       if(IC_term1(ic).eq.0) IC_term1(ic)=it
                             IC_term2(ic)=it
      End do

      End Subroutine Sort_terms

