!======================================================================
      Integer Function Ipointer (n,ii,i)
!======================================================================
!     Determines position of an element i in the integer array ii
!     A zero value indicates that the element is absent
!----------------------------------------------------------------------
      Integer, intent(in) :: n,ii(n),i

      ipointer = 0
      Do j = 1,n
        if (ii(j).eq.i) then; ipointer=j; Return; end if
      End do

      End Function Ipointer


!======================================================================
      Integer Function Rpointer (n,rr,r)
!======================================================================
!     Determines position of an element i in the set ii
!     A zero value indicates that element is not a member of the set
!----------------------------------------------------------------------
      Integer, intent(in) :: n
      Real(8), intent(in) :: r, rr(n)

      Rpointer = 0
      Do j = 1,n
        if (rr(j).eq.r) then; Rpointer=j; Return; end if
      End do

      End Function Rpointer


!======================================================================
      Integer Function Apointer (n,AA,a)
!======================================================================
!     Determines position of an element i in the set ii
!     A zero value indicates that element is not a member of the set
!----------------------------------------------------------------------
      Integer, intent(in)      :: n
      Character(*), intent(in) :: a,AA(n)

      Apointer = 0
      Do j = 1,n
        if (AA(j).eq.a) then; Apointer=j; Return; end if
      End do

      End Function Apointer
