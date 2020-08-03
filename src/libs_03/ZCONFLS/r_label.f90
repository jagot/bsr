!======================================================================
      Subroutine R_label(nu,ii,kset)
!======================================================================
!     generates LABEL list from c-file (unit nu)
!----------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer, intent(in) :: nu,ii,kset
      Integer :: nc
      
      if(allocated(Label)) Deallocate(Label)
      Allocate(Label(ncfg))

      nc=0
      rewind(nu)
    1 read(nu,'(a)',end=2) CONFIG
      if(CONFIG(1:1).eq.'*') go to 2
      if(CONFIG(5:5).ne.'(') go to 1
      read(nu,'(a)') COUPLE
      nc=nc+1
      if(nc.gt.ncfg) Stop 'R_label: nc  > ncfg'
      Call DECODE_c
      Call LABEL_C(LABEL(nc),ii,kset)
      go to 1
    2 Continue
      if(nc.ne.ncfg) Stop ' R_label: nc <> ncfg'

      End Subroutine R_label
