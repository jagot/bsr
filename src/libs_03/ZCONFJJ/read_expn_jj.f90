!====================================================================
      Subroutine read_expn_jj(nu)
!====================================================================
!     reads the expansion coefficients from GRASP c-file (unit in)
!--------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer, intent(in) :: nu
      Integer :: ic

      if(allocated(WC)) Deallocate(WC); Allocate(WC(ncfg))
      WC = 0.d0
      rewind(nu)
      ic=0
    1 read(nu,'(a)',end=2) AS
      if(AS(6:6).ne.'(') go to 1
      ic=ic+1
      if(ic.gt.ncfg) then 
        write(*,*) 'ic = ',ic
        write(*,*) 'ncfg = ',ncfg
        Stop ' read_expn_jj:  ic > ncfg'
      end if
      ia=INDEX(AS,')',BACK=.TRUE.)+1
      if(LEN_TRIM(AS(ia:)).ne.0) read(AS(ia:),*) WC(ic)
      go to 1
    2 if(ic.ne.ncfg) then
        write(*,*) 'ic = ',ic
        write(*,*) 'ncfg = ',ncfg
        Stop ' read_expn_jj: ic <> ncfg'
      end if
      End Subroutine read_expn_jj
