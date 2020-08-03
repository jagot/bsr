!======================================================================
      Subroutine ephase (n,m,K,uk,ui,us)
!======================================================================
!     n - number of open channels
!     m - number of channels (max.dimension)
!     K - symmetrical matrix
!     uk - eigenvectors of K-matrix 
!     ui - eigenphases (on module of PI) ->  tan^-1 (eigenvalue) / PI
!     us - sum of eigenphases (on module of PI)
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(in)  :: n,m
      Real(8), intent(in)  :: K(m,m)
      Real(8), intent(out) :: uk(m,m),ui(m),us
      Real(8) :: pi
      Integer :: i,info

! ... get the eigenvalues and eigenvectors of the K-matrix:

      uk=K;   Call LAP_DSYEV('V','L',n,m,uk,ui,info)
      if(info.ne.0) then
       ui=0; us=0; Return 
      end if

! ... eigenphases:

      pi = dacos(-1.d0)
      us = 0.d0
      Do i=1,n
       ui(i)  = datan(ui(i))/pi
       us = us + ui(i)
      End do

      End Subroutine ephase 

