!====================================================================
      Subroutine R_term(nu)
!====================================================================
!     reads the terms list from c-file (unit nu)
!--------------------------------------------------------------------
      Use conf_LS
      Use term_LS

      Implicit none
      Integer, intent(in) :: nu
      Integer :: nc
      Integer, external :: Idef_ncfg,Ifind_term

      Call Alloc_term_LS(0)
      Call Alloc_term_LS(iterms)
      ncfgt = Idef_ncfg(nu)
      Call Alloc_LSP(ncfgt)

      rewind(nu)
      nc = 0
    1 Read(nu,'(a)',end=2) CONFIG
      if(CONFIG(1:1).eq.'*') go to 2
      if(CONFIG(5:5).ne.'(') go to 1
      Read(nu,'(a)') COUPLE
      Call DECODE_c
      nc = nc + 1
      LSP(nc)= Ifind_term(Ltotal,Stotal,2)
      go to 1
    2 Rewind(nu)

      End Subroutine R_term

