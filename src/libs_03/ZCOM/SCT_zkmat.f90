!=======================================================================
      Subroutine ZKMAT (nch,nopen,RA,RB,RMAT,F,G,FP,GP,KMAT)
!=======================================================================
!     CALCULATE THE K-MATRIX FROM THE R-MATRIX AND ASYMPTOTIC
!     SOLUTIONS ON THE BOUNDARY:
!
!     B K = A,    A = -S + R*(a*S' - b*S),     B = C - R*(a*C' - b*C),
!
!     where S,C - so-called sin(cos) solutions and their derivatives
!                 [F,G, FP,GP]
!
!     OUTPUT: KMAT   -    K-MATRIX
!-----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nch,nopen
      Real(8), intent(in) :: RA,RB
      Real(8), intent(in) :: F(nch,nch),G(nch,nch),FP(nch,nch),GP(nch,nch)
      Real(8), intent(in) :: RMAT(nch,nch)
      Real(8), intent(out):: KMAT(nch,nch)
      Integer :: info
      Real(8) :: A(nch,nch),B(nch,nch),C(nch,nch)

      if(nopen.le.0) Stop 'ZKMAT: nopen <= 0'
      if(nopen.gt.nch) Stop 'ZKMAT: nopen > nch'

! ... cos-solutions:

      A = RA*GP - RB*G
      C = MATMUL(RMAT,A)
      B = G - C

! ... sin-solutions:

      A = RA*FP - RB*F
      C = MATMUL(RMAT,A)
      KMAT = C - F

      Call LAP_DGESV(nch,nch,nopen,B,KMAT,info)

      End Subroutine ZKMAT
