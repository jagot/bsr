!======================================================================
      Integer Function New_index(l,kksmax,nwf,ln,kn)
!======================================================================
!     assign new set index for orbital 'l'  
!----------------------------------------------------------------------
      Use param_LS, only: ksmax

      Implicit none
      Integer, Intent(in) :: l,kksmax,nwf
      Integer, Intent(in) :: ln(*),kn(*)
      Integer :: i,k,m

      Do m=1,ksmax                    
       k=m
       Do i=1,nwf
        if(l.eq.ln(i).and.m.eq.kn(i)) then
         k=0; Exit
        end if
       End do
       if(k.eq.m) Exit
      End do
      New_index = k

      End Function New_index
