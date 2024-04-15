!====================================================================
      Subroutine Jdef_JP(nu,J,P)
!====================================================================
!     define term (J,pi) from the GRASP c-file (unit nu)
!--------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer, intent(in) :: nu
      Integer, intent(out) :: J,P

      J=0; P=0
      rewind(nu)
    1 Read(nu,'(a)') CONFIG; if(CONFIG(6:6).ne.'(') go to 1
      Read(nu,'(a)') SHELLJ
      Read(nu,'(5x,a)') INTRAJ
      Call Decode_cj
      P=SUM(ln(1:no)*iq(1:no)); P=(-1)**P
      J=Jintra(no)

      End Subroutine Jdef_JP


!====================================================================
      Subroutine Jdef_JPE(nu,J,P,E)
!====================================================================
!     define state term and energy from the GRASP c-file (unit nu)
!--------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer, intent(in) :: nu
      Integer, intent(out) :: J,P
      Real(8), intent(out) :: E

      J=0; P=0; E=0.d0
      rewind(nu)
      read(nu,'(a)') AS
      if(AS(1:4).eq.'Core') then; read(AS(16:),*) E
    1 Read(nu,'(a)') CONFIG; if(CONFIG(6:6).ne.'(') go to 1
      Read(nu,'(a)') SHELLJ
      Read(nu,'(5x,a)') INTRAJ
      Call Decode_cj
      P=SUM(ln(1:no)*iq(1:no)); P=(-1)**P
      J=Jintra(no)
      end if

      End Subroutine Jdef_JPE

