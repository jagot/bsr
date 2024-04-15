!======================================================================
      Subroutine ZSMAT (nch,nopen,KMAT,SMAT)
!======================================================================
!     Convert K-matrix to S-matrix  by solving the equation:   AS = B
!     where  A=(1-iK)  and   B=(1+iK)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nch,nopen
      Real(8), intent(in) :: KMAT(nch,nch)
      Complex(8), intent(out) :: SMAT(nopen,nopen)
      Integer :: i,j,info
      Complex(8) :: Ident(nopen,nopen), A(nopen,nopen), B(nopen,nopen)
      Complex(8) :: zero = (0.d0,0.d0), one = (1.d0,0d0), &
                    onei = (0.d0,1.d0)

      if(nopen.le.0)   Stop 'ZSMAT: nopen <= 0'
      if(nopen.gt.nch) Stop 'ZSMAT: nopen > nch'

      Ident = zero;  Do i=1,nopen; Ident(i,i)=one; End do

      Do i = 1,nopen
       Do j = 1,nopen
        A(i,j) = Ident(i,j) - onei*KMAT(i,j)
        B(i,j) = Ident(i,j) + onei*KMAT(i,j)
       End do
      End do

      Call LAP_ZGESV(nopen,nch,A,B,info)

      SMAT = B

      End Subroutine ZSMAT
