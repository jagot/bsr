!======================================================================
      Subroutine read_conf_jj(muc,kshift,job,check)
!======================================================================
!     read and add configurations to the list "conf_jj"
!     job  =  'add'     -  just add
!          =  'detect'  -  return -ic if exist and 
!          =   others   -  add if not exist 
!     check = 'check'   -  check the configurations for number of 
!                          electrons and parity 
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Character(*), intent(in) :: job, check 
      Integer, intent(in) :: muc,kshift
      Integer, external :: Iadd_cfg_jj
      Integer :: nuc,i,ic

      nuc=iabs(muc); if(muc.gt.0) rewind(nuc)
      if(check.eq.'check') then; ne=0; parity=0; end if
    1 read(nuc,'(a)',end=2) CONFIG
      if(CONFIG(1:3).eq.'*') go to 2
      if(CONFIG(6:6).ne.'(') go to 1
      read(nuc,'(a)') SHELLJ
      read(nuc,'(5x,a)') INTRAJ
      Call Decode_cj
      in = in + kshift
      if(check.eq.'check') Call Test_cj
      ic = Iadd_cfg_jj(job)
      if(ic.lt.0) Stop 'Read_conf_jj: repeated states?'
      WC(ic)=0.d0
      i=INDEX(CONFIG,')',BACK=.TRUE.)+1
      if(LEN_TRIM(CONFIG(i:)).ne.0) Read(CONFIG(i:),*) WC(ic)
      go to 1
    2 Continue

      End Subroutine read_conf_jj

!======================================================================
      Subroutine read_config_jj(nuc)
!======================================================================
!     Read only configurations from GRASP c-file (unit 'nuc')
!     (no weights)
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer, intent(in) :: nuc
      Integer, external :: Iadd_cfg_jj
      Integer :: i

      SHELLJ =' '
      INTRAJ =' '
      Jshell=0; Vshell=0; Jintra=0

      rewind(nuc); ne=0; parity=0
    1 read(nuc,'(a)',end=2) CONFIG
      if(CONFIG(1:1).eq.'*') go to 2
      if(CONFIG(6:6).ne.'(') go to 1
      Call Decode_confj
      i = Iadd_cfg_jj('find')
      go to 1
    2 Continue

      End Subroutine read_config_jj


