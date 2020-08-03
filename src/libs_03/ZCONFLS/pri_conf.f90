!======================================================================
      Subroutine Pri_conf (nu,ic,w)
!======================================================================
!     prints 'ic' configuration from module 'conf_LS'
!----------------------------------------------------------------------
      Use conf_LS

      Integer, intent(in) :: nu
      Real(8), intent(in) :: w

      if(ic.lt.0) Stop ' Pri_ic: ic < 0 '
      if(ic.gt.ncfg) Stop ' Pri_ic: ic > ncfg '

      if(ic.gt.0) Call Get_cfg_LS(ic)

      Call Incode_c
      if(w.ne.0.d0) then
       write(nu,'(a64,F11.8)') CONFIG,w
      else
       write(nu,'(a)') trim(CONFIG)
      end if
      write(nu,'(a)') trim(COUPLE)

      End Subroutine Pri_conf


!======================================================================
      Subroutine Prj_conf_LS (nu,w)
!======================================================================
!     prints configuration from module 'conf_LS'
!----------------------------------------------------------------------
      Use conf_LS

      Integer, intent(in) :: nu
      Real(8), INTENT(in) :: w

      Call Incode_c
      if(w.ne.0.d0) then
       write(nu,'(a64,F11.8)') CONFIG,w
      else
       write(nu,'(a)') trim(CONFIG)
      end if
      write(nu,'(a)') trim(COUPLE)

      End Subroutine Prj_conf_LS

