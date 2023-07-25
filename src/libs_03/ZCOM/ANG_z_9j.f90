!======================================================================
      Real(8) Function Z_9j (j1,j2,j3,j4,j5,j6,j7,j8,j9)
!======================================================================
!
!     determination of 9j-symbols     
!                                    
! {j1 j2 j3}                           {j1 j2 j3} {j4 j5 j6} {j7 j8 j9}        
! {j4 j5 j6} = SUM(j) (-1)^(2j) (2j+1) 
! {j7 j8 j9}                           {j6 j9 j } {j2 j  j8} {j  j1 j4}
!
!
!     (the momenta are used in (2J+1)-representation)
!
!----------------------------------------------------------------------

      Implicit none
      Integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7,j8,j9
      Integer :: j,i1,i2
      Real(8), External :: Z_6j
      
      i1 = max(abs(j1-j9),abs(j2-j6),abs(j4-j8))+1
	  i2 = min(j1+j9,j2+j6,j4+j8)-1

	  Z_9j = 0.d0; if(i1.gt.i2) Return;  if(mod(i2-i1,2).ne.0) Return

      Do j = i1,i2,2
       Z_9j = Z_9j + j * (-1)**(j-1) * &
                     Z_6j(j1,j2,j3,j6,j9,j ) * &
                     Z_6j(j4,j5,j6,j2,j ,j8) * &
                     Z_6j(j7,j8,j9,j ,j1,j4) 
      End do
       
      End Function Z_9j


!======================================================================
      Real(8) Function Z_9jj (j1,j2,j3,j4,j5,j6,j7,j8,j9)
!======================================================================

      Implicit none
      Integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7,j8,j9
      Real(8), external :: Z_9j

      z_9jj = Z_9j (j1+j1+1,j2+j2+1,j3+j3+1,j4+j4+1,j5+j5+1,j6+j6+1, &
                    j7+j7+1,j8+j8+1,j9+j9+1)

      End Function Z_9jj
