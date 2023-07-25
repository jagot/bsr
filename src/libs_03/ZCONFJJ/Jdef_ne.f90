!====================================================================
      Integer Function Jdef_ne(nu)
!====================================================================
!     defines the number of electrons from c-file (unit nu)
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(450) ::  LINE
      Integer :: i,ii,ip,iq,ne

      ne = 0
      rewind(nu)
    1 Read(nu,'(a)',end=2) LINE
      if(LINE(1:3).eq.'***') go to 2
      if(LINE(6:6).ne.'(') go to 1
      ii = INDEX(LINE,')',BACK=.TRUE.)/9
      Do i = 1,ii
       ip = (i-1)*9+7
       read(LINE(ip:ip+1),'(i2)') iq
       ne = ne + iq
      End do
    2 rewind(nu)
      Jdef_ne = ne

      End Function Jdef_ne
