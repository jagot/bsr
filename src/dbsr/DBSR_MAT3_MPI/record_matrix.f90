!======================================================================
      Subroutine Record_matrix
!======================================================================
!     record the overlap or Hamiltonian matrix to unit nui
!----------------------------------------------------------------------
      Use dbsr_mat
#ifdef DEBUG_SPEEDUPS
      Use Timer
#endif

      Implicit none
      Integer :: i,j,k, ni,nj, ip,jp, ich,jch
      Real(8) :: S, w(ns,ms), c(ns,ns)

#ifdef DEBUG_SPEEDUPS
      Real(8) :: w_ref(ms), c_ref(ns,ns)
      Real(8) :: discrepancy
      Real(8), parameter :: tolerance = sqrt(epsilon(1.d0))
#endif

! ... non-diagonal channel blocks:

      Do ich = 2,nch;    ni=ipsol(ich)
      Do jch = 1,ich-1;  nj=ipsol(jch)
       k=icc(ich,jch); if(k.eq.0) Cycle
       if(maxval(abs(hch(:,:,k))).lt.1.d-20) Cycle
       c = 0.d0

#ifdef DEBUG_SPEEDUPS
       call TimerStart('Record_matrix: new way 1')
#endif
       w(1:ni,:) = matmul(transpose(diag(:,1:ni,ich)), hch(:,:,k))
       c(1:ni,1:nj) = matmul(w(1:ni,:), diag(:,1:nj,jch))

#ifdef DEBUG_SPEEDUPS
       call TimerStop('Record_matrix: new way 1')

       c_ref = 0.d0
       call TimerStart('Record_matrix: old way 1')
       Do i=1,ni
          Do j=1,ms
             w_ref(j)  = SUM(diag(:,i,ich)*hch(:,j,k))
          End do
          Do j=1,nj
             c_ref(i,j)= SUM(w_ref(:)*diag(:,j,jch))
          End do
       End do
       call TimerStop('Record_matrix: old way 1')

       discrepancy = sum(abs(c_ref(1:ni,1:nj)-c(1:ni,1:nj)))
       write(*,'("Record_matrix, discrepancy 1: ",e26.16," ",e26.16)') &
            discrepancy
       if(discrepancy > tolerance) error stop "Discrepancy exceeds tolerance"
#endif

       write(nui) ich,jch
       write(nui) c(1:ni,1:nj)
       if(icase.eq.0) overlaps(ich,jch) = maxval(abs(c(1:ni,1:nj)))
      End do; End do

! ... pertubers:

      if(npert.gt.0) then

! ... channel-perturber rows:

      Do ich = 1,nch;    ni = ipsol(ich)
      Do ip  = 1,npert
       k = icb(ich,ip); if(k.eq.0) Cycle
       S = maxval(abs(hcp(:,k)))
       if(S.lt.1.d-20) Cycle
       w = 0.d0

#ifdef DEBUG_SPEEDUPS
       call TimerStart('Record_matrix: new way 2')
#endif
       w(1:ni,1) = matmul(transpose(diag(:,1:ni,ich)), hcp(:,k))

#ifdef DEBUG_SPEEDUPS
       call TimerStop('Record_matrix: new way 2')
       w_ref = 0.d0

       call TimerStart('Record_matrix: old way 2')
       Do i=1,ni
          w_ref(i)=SUM(diag(:,i,ich)*hcp(:,k))
       End do
       call TimerStop('Record_matrix: old way 2')

       discrepancy = sum(abs(w_ref(1:ni)-w(1:ni,1)))
       write(*,'("Record_matrix, discrepancy 2: ",e26.16)') &
            discrepancy
       if(discrepancy > tolerance) error stop "Discrepancy exceeds tolerance"
#endif

       if(icase.eq.0) overlaps(nch+ip,ich) = maxval(abs(w(1:ni,1)))
       write(nui) nch+ip,ich
       write(nui) w(1:ni,1)
      End do; End do

! ... perturter-perturber elements:

      Do ip = 1,npert; Do jp = 1,ip; k = ibb(ip,jp)
       if(k.eq.0) Cycle
       if(hp(k).eq.0.d0) Cycle
       write(nui) ip+nch,jp+nch
       write(nui) hp(k)
       if(icase.eq.0) overlaps(nch+ip,nch+jp) = abs(hp(k))
      End do; End do

      end if   ! over npert

      write(nui) 0,0  ! sign of the end

      End Subroutine Record_matrix
