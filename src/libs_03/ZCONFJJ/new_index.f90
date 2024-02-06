!======================================================================
      Integer Function new_index(l,ksmax,nwf,ln,kn)
!======================================================================
!     assign a set index for orbital 'nl' to be different from already
!     existing ones in the list {ln,kn)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l,ksmax,nwf
      Integer, intent(in) :: ln(nwf),kn(nwf)
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

      End Function new_index
