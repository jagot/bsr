!====================================================================
      Subroutine R_CLOSED (nu)
!====================================================================
!     reads the common cloused shells list from c-file (unit nu)
!--------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer, intent(in) :: nu
      Character(4) :: EL
      Integer, external :: Ifind_nlk
      Integer :: i,j, n,l,k

      rewind(nu)
      read(nu,*)
      read(nu,'(a)') CLOSED

      NCLOSD=0
      i=1
      Do
       EL=CLOSED(i:i+3)
       if(LEN_TRIM(EL).eq.0) Exit
       NCLOSD=NCLOSD+1
       Call EL4_nlk(EL,n,l,k)
       j=Ifind_nlk(n,l,k,2)
       i=i+4
      End do
      ncore = nclosd

      End Subroutine R_CLOSED
