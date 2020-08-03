!======================================================================
      Subroutine Prepare_it_oper(eps_soo)
!----------------------------------------------------------------------
      Use param_LS
      Use symc_list_LS
      Use symt_list_LS
      Use conf_LS

      Implicit none
      Real(8) :: eps_soo
      Real(8), allocatable :: C_term(:)
      Integer :: i,j, ik,jk, ic,jc, it,jt
      Integer(8) :: ijc
      Integer(8), external :: Def_ij8

!-------------------------------------------------------------------------		 
! ... additional restrictions for spi_orbit interaction:

      if(eps_soo.gt.0.d0) then
       Allocate(C_term(nsymt));  C_term = 0.d0
       Do i=1,nsymt; ik= IT_state1(i); jk=IT_state2(i) 
        Do  j = ik,jk; ic=IP_stat(j)
          C_term(i) = max(C_term(i),abs(WC(ic)))
        End do
       End do
      end if

!-------------------------------------------------------------------------		 
! ... allocate IT_oper:

      Call Alloc_it_oper_LS(1)

!-------------------------------------------------------------------------		 
! ... define strong orthogonality:      

      Do ic = 1,nsymc
       Call Get_symc_LS(ic,Ltotal1,Stotal1,no1,nn1,ln1,iq1,kn1) 
      Do jc = 1,ic
       Call Get_symc_LS(jc,Ltotal2,Stotal2,no2,nn2,ln2,iq2,kn2) 
      
       Call Jort_conf(joper)

       Do ik=IC_term1(ic),IC_term2(ic);  it=IT_sort(ik)
       Do jk=IC_term1(jc),IC_term2(jc);  jt=IT_sort(jk)

         Do i=1,noper
          if(ioper(i).eq.0) then
           koper(i)=-1
          else
           koper(i) = joper(i)            
          end if
         End do

         if(eps_soo.gt.0.d0) then
          if(C_term(it)*C_term(jt).lt.eps_soo) koper(4:7)=-1
         end if

         ij = DEF_ij8(it,jt)
         IT_oper(:,ij) = koper(:)
       End do
       End do

      End do
      End do

      if(allocated(C_term)) Deallocate(C_term)

!-------------------------------------------------------------------------		 
! ... define IT_need and JT_need:

      if(allocated(IT_need)) Deallocate(IT_need)
                             Allocate(IT_need(nsymt))
      if(allocated(JT_need)) Deallocate(JT_need)
                             Allocate(JT_need(ij_oper))

       IT_need=0; JT_need=0
       Do it=1,nsymt; Do jt=1,it; ij=DEF_ij8(it,jt)
       Do i=1,noper
         if(IT_oper(i,ij).eq.1) Cycle
         if(IT_oper(i,ij).eq.-1) Cycle
         JT_need(ij)=1; IT_need(it)=1; IT_need(jt)=1
        End do
       End do; End do

!-------------------------------------------------------------------------		 
! ... define IC_need and JC_need:

      if(allocated(IC_need)) Deallocate(IC_need)
                             Allocate(IC_need(nsymc))
      if(allocated(JC_need)) Deallocate(JC_need)
                             ijc = DEF_ij8(nsymc,nsymc)
                             Allocate(JC_need(ijc))

      IC_need = 0; JC_need = 0
      Do ic = 1,nsymc; Do jc = 1,ic; ijc=DEF_ij8(ic,jc)
       Do ik=IC_term1(ic),IC_term2(ic);  it=IT_sort(ik)
       Do jk=IC_term1(jc),IC_term2(jc);  jt=IT_sort(jk)
        ij = DEF_ij8(it,jt)
        if(JT_need(ij).eq.0) Cycle
        JC_need(ijc)=1; IC_need(ic)=1; IC_need(jc)=1
       End do
       End do
      End do; End do

      Deallocate(JT_need)

      End Subroutine Prepare_it_oper



!======================================================================
      Subroutine Record_it_oper(nu)
!======================================================================
      Use symc_list_LS, only: nsymc,IC_need,JC_need
      Use symt_list_LS, only: nsymt,ij_oper,IT_oper,IT_need

      Implicit none
      Integer, intent(in) :: nu

      write(nu) ij_oper,nsymt,nsymc
      write(nu) IT_oper
      write(nu) IT_need
      write(nu) IC_need
      write(nu) JC_need
            
      End Subroutine Record_it_oper


!======================================================================
      Subroutine Load_it_oper(nu)
!======================================================================
      Use symc_list_LS, only: nsymc,IC_need,JC_need
      Use symt_list_LS, only: noper,nsymt,ij_oper,IT_oper,IT_need

      Implicit none
      Integer, intent(in) :: nu
      Integer :: it,ic
      Integer(8) :: ijc
      Integer(8), external :: Def_ij8

      read(nu) ij_oper,it,ic
      if(ic.ne.nsymc) Stop 'Load_it_oper:  nsymc - ?'
      if(it.ne.nsymt) Stop 'Load_it_oper:  nsymt - ?'

      if(allocated(IT_oper)) Deallocate(IT_oper)
      Allocate(IT_oper(noper,ij_oper))
      read(nu) IT_oper

      if(allocated(IT_need)) Deallocate(IT_need)
      Allocate(IT_need(nsymt))
      read(nu) IT_need

      if(allocated(IC_need)) Deallocate(IC_need)
      Allocate(IC_need(nsymc))
      read(nu) IC_need

      if(allocated(JC_need)) Deallocate(JC_need)
      ijc = DEF_ij8(nsymc,nsymc)
      Allocate(JC_need(ijc))
      read(nu) JC_need
            
      End Subroutine Load_it_oper

