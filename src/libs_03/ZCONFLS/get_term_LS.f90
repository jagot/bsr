!======================================================================
      Subroutine Term_ic (ic,L,S)
!======================================================================
!     term of the state ic
!----------------------------------------------------------------------

      Use conf_LS
      Use symc_list_LS
      Use symt_list_LS

      Implicit none
      Integer :: ic,L,S

      iterm=IC_term(ic)
      iconf=it_conf(iterm)
      L = LT_conf(iconf)
      S = ST_conf(iconf)

      End Subroutine Term_ic
