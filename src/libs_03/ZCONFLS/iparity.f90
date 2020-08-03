!====================================================================
      Integer Function Iparity()
!====================================================================
!     defines parity for configuration in module 'conf_LS'
!--------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer :: i,m

      m=0; Do i=1,no;  m=m+ln(i)*iq(i);  End do
      Iparity=(-1)**m

      End Function Iparity

!====================================================================
      Subroutine Conf_parity(k)
!====================================================================
!     defines parity for ZOI configuration in module 'conf_LS'
!--------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer :: k,i,m

      m = 0
      if(k.eq.0) then
       Do i=1,no;  m=m+ln(i)*iq(i);  End do; Ptotal=(-1)**m
      elseif(k.eq.1) then
       Do i=1,no1;  m=m+ln1(i)*iq1(i);  End do; Ptotal1=(-1)**m
      elseif(k.eq.2) then
       Do i=1,no2;  m=m+ln2(i)*iq2(i);  End do; Ptotal2=(-1)**m
      end if

      End Subroutine Conf_parity


!====================================================================
      Integer Function Def_ne()
!====================================================================
!     define number of electrons in configuration in module 'configs'
!--------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer :: i,m

      m=0; Do i=1,no;  m=m+iq(i);  End do
      Def_ne=m

      End Function Def_ne

