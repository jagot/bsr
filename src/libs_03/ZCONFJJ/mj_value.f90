!====================================================================
      Integer Function Index_mj (mj)
!====================================================================
!     index of jm-orbital in the common list: -1,+1,-3,+3,-5,+5, ...
!     mj -> in 2j-representation  (odd integers)
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: mj
      Index_mj = iabs(mj)
      if(mj.gt.0) index_mj = index_mj + 1
      END Function Index_mj

!====================================================================
      Integer Function mj_value(i)
!====================================================================
!     mj value for orbital 'i' in the lit: -1,+1,-3,+3,-5,+5, ...
!     mj -> in 2j-representation
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: i
      mj_value = -i
      if(mod(i,2).eq.0) mj_value = i - 1
      END Function mj_value

!====================================================================
      Integer Function ndets_jq(j,q)
!====================================================================
!     number of det.s in subshells j^k (Newton's binom)
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j,q
      Integer :: i
      Real(8) :: S
      if(q.gt.j+1) Stop 'ndets_jq:  q > q_max'
      S=1.d0
      Do i=q+1,j+1; S=S*i/(i-q); End do
      ndets_jq = S + 0.1d0
      End Function ndets_jq

