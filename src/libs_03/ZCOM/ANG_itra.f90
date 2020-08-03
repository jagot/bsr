!====================================================================
      Integer Function ITRA (i1,i2,i3)
!====================================================================
!     check triangle relations for  momentums i1,i2,i3
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: i1,i2,i3

      ITRA=1
      if(i1.gt.i2+i3.or.i2.gt.i1+i3.or.i3.gt.i1+i2) ITRA=0
      if(i1.lt.iabs(i2-i3).or.i2.lt.iabs(i1-i3).or. &
         i3.lt.iabs(i1-i2)) ITRA=0

      End function ITRA


!====================================================================
      Integer Function ITRI (i1,i2,i3)
!====================================================================
!     check triangle relations for the i1,i2,i3 momentums given
!     in the (2J+1)-format
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: i1,i2,i3

      ITRI=1
      if(i1.gt.i2+i3-1.or.i2.gt.i1+i3-1.or.i3.gt.i1+i2-1) ITRI=0
      if(i1.lt.iabs(i2-i3)+1.or.i2.lt.iabs(i1-i3)+1.or. &
         i3.lt.iabs(i1-i2)+1) ITRI=0

      End Function ITRI

                                       