!======================================================================
      Subroutine Record_matrix(nu)
!======================================================================
!     record interaction/overlap matrix to the unit 'nu'
!----------------------------------------------------------------------
      Use bsr_mat

      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,ip,jp,ich,jch
      Real(8) :: S

! ... channel-channal blocks:

      Do ich=1,nch; Do jch=1,ich;  i=icc(ich,jch);  if(i.eq.0) Cycle
       S = SUM(abs(hcc(:,:,i)))
       if(S.eq.0.d0) Cycle
       write(nu) ich,jch
       write(nu) hcc(:,:,i)
      End do; End do

! ... pertubers:

      if(npert.gt.0) then

       Do ich=1,nch; Do ip=1,npert; i=icb(ich,ip); if(i.eq.0) Cycle
        S = SUM(abs(hcb(:,i)))
        if(S.eq.0.d0) Cycle
        write(nu) ip+nch,ich
        write(nu) hcb(:,i)
       End do; End do

       Do ip=1,npert; Do jp=1,ip; i=ibb(ip,jp); if(i.eq.0) Cycle
        if(hbb(i).eq.0.d0) Cycle
        write(nu) ip+nch,jp+nch
        write(nu) hbb(i)
       End do; End do

      end if

! ... sign of the end

      write(nu) -1,-1

      End Subroutine Record_matrix

