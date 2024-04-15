!====================================================================
      Integer Function jopen (e,ilsp)
!====================================================================
!     number of open channels in the jk-case
!--------------------------------------------------------------------
      Use target
      Use channels

      Implicit  none
      Real(8) :: e
      Integer :: ilsp,i

      jopen = 0
      Do i = 1,nch(ilsp)
       if(e.lt.etarg(iptar(ilsp,i))) Exit
       jopen = jopen + 1
      End do

      End Function jopen


!====================================================================
      Integer Function jopen_ion (e,ilsp)
!====================================================================
!     number of open channels
!--------------------------------------------------------------------
      Use target_ion
      Use channels_ion

      Implicit  none

      Real(8) :: e
      Integer :: ilsp,i

      jopen_ion = 0
      Do i = 1,nch(ilsp)
       if(e.lt.etarg(iptar(ilsp,i))) Exit
       jopen_ion = jopen_ion + 1
      End do

      End Function jopen_ion
