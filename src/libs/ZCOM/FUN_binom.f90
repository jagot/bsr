!====================================================================
      Integer Function Iglq(l,iq)
!====================================================================
!     total stat.weight of lq-subshell (Newton's binom)
!--------------------------------------------------------------------
      kl=4*l+2
      if(iq.gt.kl) Stop 'Iglq:  iq > max.'
      S=1.0
      Do i=iq+1,kl
       S=S*i/(i-iq)
      End do
      Iglq=S+0.5

      End Function Iglq

