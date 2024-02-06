!======================================================================
      Subroutine Check_conf_jj(nuc,nub)
!======================================================================
!     Read the configuration list from c-file (unit 'nuc'),
!     define there angular symmetries and compare with ones in the
!     angular coefficients data bank (unit 'nub')
!     Prepare other angular arrays.
!----------------------------------------------------------------------
      Use conf_jj
      Use symc_list
      Use symt_list

      Implicit none
      Integer, intent(in) :: nuc,nub
      Integer, external :: Jdef_ncfg, Iadd_cfg_jj
      Integer :: i,k,it,jt,ij,nsymc0,nsymt0

! ... define ncfg:

      Call Alloc_cfg(0)
      Call Alloc_symc(0)
      Call Alloc_symt(0)

!      ncfg=Jdef_ncfg(nuc);
!      if(ncfg.eq.0) Stop ' ncfg = 0, nothing to do '

! ... read bank information:

      Call Read_symc(nub); nsymc0 = nsymc
      Call Read_symt(nub); nsymt0 = nsymt

!---------------------------------------------------------------------
! ... define symmetries from c-file:

      rewind(nuc); ne=0; parity=0
    1 read(nuc,'(a)',end=2) CONFIG
      if(CONFIG(1:3).eq.'***') go to 2
      if(CONFIG(6:6).ne.'(') go to 1
      read(nuc,'(a)') SHELLJ
      read(nuc,'(5x,a)') INTRAJ
      Call Decode_cj;  Call Test_cj
      i = Iadd_cfg_jj('add')
      go to 1
    2 Continue

      if(nsymt.gt.nsymt0.or.nsymc.gt.nsymc0) &
       Stop 'bnk-file is not complete; run DBSR_BREIT first!'

!----------------------------------------------------------------------
! ... define IT_state:

      if(allocated(IS_order )) Deallocate(IS_order)
      if(allocated(IT_state1)) Deallocate(IT_state1)
      if(allocated(IT_state2)) Deallocate(IT_state2)
      Allocate(IS_order(ncfg),IT_state1(nsymt),IT_state2(nsymt))
      IT_state1=0; IT_state2=0

      Call SORTI (ncfg,IS_term,IS_order)

      Do i=1,ncfg
       it=IS_term(IS_order(i))
       if(IT_state1(it).eq.0) IT_state1(it)=i; IT_state2(it)=i
      End do

!----------------------------------------------------------------------
! ... define if we need additional angular coefficients:

      if(allocated(IT_done)) Deallocate(IT_done)
                             Allocate(IT_done(nsymt*(nsymt+1)/2))
      Call Read_done(nub)

      k = 0
      Do it = 1,nsymt;  if(IT_state1(it).eq.0) Cycle
       Do jt = 1,it;    if(IT_state1(jt).eq.0) Cycle
        ij = (it-1)*it/2+jt
        if(IT_done(ij).gt.0) Cycle
        k=1; Exit
       End do
        if(k.eq.1) Exit
      End do
      if(k.eq.1) Stop 'int_bnk file is not complete; run DBSR_BREIT first'

      Deallocate(IT_done)

      End Subroutine Check_conf_jj


