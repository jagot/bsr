!======================================================================
      Subroutine Check_file(AF)
!======================================================================
!     stop with message if given file AF does not exist
!----------------------------------------------------------------------
      Character(*), Intent(in) :: AF
      Logical :: EX
      Character(256) :: error_msg
      Inquire(file=trim(AF),exist=EX)
      if(.not.EX) then
         write(error_msg, '(" can not find file  ",a)') trim(AF)
         Error Stop error_msg
      end if
      End Subroutine Check_file


!======================================================================
      Integer Function Icheck_file(AF)
!======================================================================
!     check if the file AF exists
!----------------------------------------------------------------------
      Character(*), Intent(in) :: AF
      Logical :: EX
      Inquire(file=trim(AF),exist=EX)
      Icheck_file = 1
      if(.not.EX) Icheck_file = 0
      End Function Icheck_file


!======================================================================
      Subroutine Find_free_unit(nu)
!======================================================================
!     provide free unit to open new file
!----------------------------------------------------------------------
      Implicit none
      Integer :: nu,i
      Logical :: connected
      nu = 0
      Do i=21,999
       INQUIRE(UNIT=i,OPENED=connected)
       if(connected) Cycle
       nu = i
       Exit
      End do
      if(nu.eq.0) Error Stop 'Find_free_unit: nu = 0'
      End Subroutine Find_free_unit
