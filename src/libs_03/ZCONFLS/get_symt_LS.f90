!======================================================================
      Subroutine Get_symt_LS(it,iconf,no,LS)
!======================================================================
!     extrats term symetry 'it'
!----------------------------------------------------------------------

      Use param_LS, only: msh
      Use symt_list_LS
      Use symc_list_LS

      Implicit none
      Integer :: it,iconf, no, i,j,ip
      Integer :: LS(msh,5)

      if(it.le.0.or.it.gt.nsymt) Stop 'Get_symt_LS: <it> is out of range'

      iconf = IT_conf(it)
      no = no_conf(iconf)

      ip = ip_term(it)
      Do i=1,no; ip=ip+1;  LS(i,:) = LS_term(ip,:);  End do

      End Subroutine Get_symt_LS