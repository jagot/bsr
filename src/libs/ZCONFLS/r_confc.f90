!====================================================================
      Subroutine R_conf_LS(nu,kshift)
!====================================================================
!     reads configurations from c-file (unit in) and saves them in
!     module conf_LS
!     the dublicated configurations are omitted
!--------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer, intent(in) :: nu
      Integer, intent(in), optional :: kshift
      Integer :: ic
      Integer, external :: Ifind_cfg_LS
      Character(100) :: AS
      Real(8) :: W

      if(mcfg.eq.0) Call alloc_cfg_LS(icfg)

      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      read(AS,'(a)') CONFIG
      W = 0.d0
      if(LEN_TRIM(AS).gt.64) read(AS(65:),*) W
      read(nu,'(a)') COUPLE
      Call Decode_c
      kn=kn+kshift
      ic = Ifind_cfg_LS()
      WC(ic) = W
      go to 1
    2 Continue

      End Subroutine R_conf_LS


!====================================================================
      Subroutine Add_conf_LS(mu,kshift)
!====================================================================
!     reads configurations from c-file (unit in) and saves them in
!     module conf_LS (with dublicated configurations are also saved)
!     if mu > 0  -  read from the beginning of the file
!     if mu < 0  -  read from the given position
!--------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer, intent(in) :: mu
      Integer, intent(in), optional :: kshift
      Integer :: ic,nu
      Integer, External :: Iadd_cfg_LS
      Character(100) :: AS
      Real(8) :: W

      if(mcfg.eq.0) Call alloc_cfg_LS(icfg)

      parity = 0; ne=0
      nu = iabs(mu); if(mu.gt.0) rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      read(AS,'(a)') CONFIG
      W = 0.d0
      if(LEN_TRIM(AS).gt.64) read(AS(65:),*) W
      read(nu,'(a)') COUPLE
      Call Decode_c
      kn=kn+kshift
      ic = Iadd_cfg_LS()
      WC(ic) = W
      Call Test_c
      go to 1
    2 Continue

      End Subroutine Add_conf_LS


!====================================================================
      Subroutine Load_conf_LS(mu,kshift)
!====================================================================
!     reads configurations from c-file (unit in) and saves them in
!     module conf_LS (with dublicated configurations are also saved)
!     if mu > 0  -  read from the beginning of the file
!     if mu < 0  -  read from the given position
!--------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer, intent(in) :: mu
      Integer, intent(in), optional :: kshift
      Integer :: ic,nu
      Integer, external :: Iadd_cfg_LS
      Character(100) :: AS
      Real(8) :: W

      if(mcfg.eq.0) Call alloc_cfg_LS(icfg)

      parity = 0; ne=0
      nu = iabs(mu); if(mu.gt.0) rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      read(AS,'(a)') CONFIG
      W = 0.d0
      if(LEN_TRIM(AS).gt.64) read(AS(65:),*) W
      read(nu,'(a)') COUPLE
      Call Decode_c
      kn=kn+kshift
      ic = Iadd_cfg_LS()
      WC(ic) = W
      go to 1
    2 Continue

      End Subroutine Load_conf_LS


!====================================================================
      Subroutine R_config_LS(nu,kshift)
!====================================================================
!     reads only configurations (without terms)
!--------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer, intent(in) :: nu
      Integer, intent(in), optional :: kshift
      Integer :: i,ic
      Integer, External :: Ifind_cfg_LS
      Character(100) :: AS
      Real(8) :: W

      if(mcfg.le.0) Call alloc_cfg_LS(icfg)

      Do i=1,msh; COUPLE((i-1)*4+1:i*4)=' 1S '; End do

      rewind(nu)
      AS = ' '
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      read(AS,'(a)') CONFIG
      W = 0.d0
      if(LEN_TRIM(AS).gt.64) read(AS(65:),*) W
      Call Decode_c
      kn=kn+kshift
      if(kshift.eq.-1) kn=0
      ic = Ifind_cfg_LS()
      if(abs(W).gt.WC(ic)) WC(ic) = abs(W)
      go to 1
    2 Continue

      End Subroutine R_config_LS







