!======================================================================
      Real(8) Function F_gaus(x,y,s)
!======================================================================
!     gives value of Gaussian with width 's' and center 'y'
!     at the point 'x'
!----------------------------------------------------------------------
      Implicit none
      Real(8), intent(in) :: x,y,s
      Real(8), save :: A,B,PI
      Integer, save :: start = 0

      if(start.eq.0) then
       PI = ACOS(-1.0)
       A  = 2.d0*sqrt(LOG(2.d0)/PI)
       B  = -4.d0*LOG(2.d0)
      end if

      F_gaus = A/S * EXP( B * ((x-y)/s)**2 )

      End Function F_Gaus
