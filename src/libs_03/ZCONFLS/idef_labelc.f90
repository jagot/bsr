!======================================================================
      Subroutine Idef_LabelC(nu,ii,kset,LabelC)
!======================================================================
!     generates LABEL list from c-file (unit nu)
!----------------------------------------------------------------------
      Use conf_LS
      Implicit none
      Integer, intent(in) :: nu,ii,kset
      Character(64), intent(out) :: LabelC
      Real(8) :: C, CC
      Integer :: i 
      Character(275) :: A
      
      LabelC = ' '

      rewind(nu)
      read(nu,'(a)',end=2) A
      i=INDEX(A,':')
      if(i.ne.0) read(A(i+2:),'(a)') LabelC
      if(len_trim(LabelC).ne.0) Return
	  
      C = 0.d0
    1 read(nu,'(a)',end=2) A
      if(A(1:1).eq.'*') go to 2
      if(A(5:5).ne.'(') go to 1
      read(A,'(a64,f11.8)') CONFIG, CC
      read(nu,'(a)') COUPLE
      if(abs(CC).lt.C) go to 1
      C = abs(CC)
      Call DECODE_c
      Call LABEL_C(LabelC,ii,kset)
      go to 1
    2 Continue


      End Subroutine Idef_LabelC
