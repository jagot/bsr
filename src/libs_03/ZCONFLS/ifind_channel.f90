!======================================================================
      Integer Function Ifind_channels_dc(ic,klsp)
!======================================================================
!     find channel (or perturber) for given configuration 'ic'
!----------------------------------------------------------------------
      Use target_dc; Use channels_dc

      Implicit none
      Integer, intent(in) :: ic,klsp
      Integer :: i,ich,kch,ip

      i=1
      ip = ipch(klsp); kch = nch(klsp)
      Do ich = 1,kch
       if(ic.gt.ipconf(ip+ich)) i=ich+1
      End do
      if(ic.gt.ipconf(ip+kch)) i=ic-ipconf(ip+kch)+kch
      Ifind_channels_dc = i

      End Function Ifind_channels_dc
