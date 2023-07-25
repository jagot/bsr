!======================================================================
      Subroutine Fano_tma (l1,l2,ar,ai,TR,TI)
!======================================================================
!     (tr,ti) =  i^(l1-l2) * (ar,ai)
!
!     used for FANO-phase correction of matrix elements
!----------------------------------------------------------------------
      Implicit none
      Real(8), intent(in)  :: ar,ai
      Real(8), intent(out) :: TR,TI
      Integer, intent(in)  :: l1,l2
      Integer :: ll, k
      
      ll = l1-l2; k = mod(ll,4)

      Select case(k)
       Case( 0); TR =  ar; TI =  ai
       Case( 1); TR = -ai; TI =  ar
       Case(-1); TR =  ai; TI = -ar
       Case( 2); TR = -ar; TI = -ai
       Case(-2); TR = -ar; TI = -ai
       Case( 3); TR =  ai; TI = -ar
       Case(-3); TR = -ai; TI =  ar
      End Select

      End Subroutine Fano_tma



 
