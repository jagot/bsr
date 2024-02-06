!====================================================================
      Subroutine print_conf_jj (nu,ic,w)
!====================================================================
!     prints the  CAS "ic" stored in the module conf_jj
!     to unit 'nu' in spectroscopic format
!     w - optional expansion coefficient
!--------------------------------------------------------------------
      Use conf_jj

      Integer, intent(in) :: nu
      Real(8), intent(in) :: w

      if(ic.gt.ncfg) Stop ' Pri_ic: ic > ncfg '

      if(ic.gt.0)  Call Get_cfg_jj(ic)

      Call Incode_cj

      ia = len_trim(CONFIG)
      if(w.ne.0.d0) then
       if(ia.le.72) then
        write(nu,'(a72,F20.15)') CONFIG(1:72),w
       else
        write(nu,'(a,F20.15)') CONFIG(1:ia),w
       end if
      else
       write(nu,'(a)') CONFIG(1:ia)
      end if

      write(nu,'(a)') SHELLJ(1:ia)
      write(nu,'(9x,a)') INTRAJ(1:ia)

      End Subroutine print_conf_jj


