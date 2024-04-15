!====================================================================
      Subroutine R_expn(in,ncfg,C)
!====================================================================
!     reads the expansion coefficients from c-file (unit in)
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: in, ncfg
      Real(8), intent(out) :: C(ncfg)
      Character  A*100
      Integer :: ic

      rewind(in)
      ic=0
    1 read(in,'(a)',end=2) A
      if(A(5:5).ne.'(') go to 1
      ic=ic+1
      if(ic.gt.ncfg) Stop ' R_expn:  ic > ncfg'
      read(A(65:),*) C(ic)
      go to 1
    2 if(ic.ne.ncfg) Stop ' R_expn: ic <> ncfg'

      End Subroutine R_expn
