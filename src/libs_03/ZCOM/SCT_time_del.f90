!======================================================================
      Subroutine time_del (nch,nopen,AK1,AK2,de,tdvalue,parg)
!======================================================================
!     time delay from two K-matrices
!----------------------------------------------------------------------
      Implicit double precision (A-H,O-Z)
      Integer, intent(in) :: nch,nopen
      Real(8), intent(in) :: AK1(nch,nch),AK2(nch,nch)
      Real(8), intent(in) :: de
      Real(8), intent(out) :: tdvalue
      Real(8), intent(out) :: parg(nopen)
      Complex(8), allocatable :: SMAT1(:,:),SMAT2(:,:),SMAT(:,:),DES(:,:),Q(:,:)
      Real(8), allocatable :: qval(:)
      Integer :: info

      Allocate( SMAT1(nopen,nopen), SMAT2(nopen,nopen), &
                SMAT(nopen,nopen), DES(nopen,nopen),  &
                Q(nopen,nopen), qval(nopen))

      Call ZSMAT(nch,nopen,AK1,SMAT1)
      Call ZSMAT(nch,nopen,AK2,SMAT2)

      SMAT = (smat1 + smat2) / 2.d0
      DES  = (smat2 - smat1) / de
      Q = (0,-1) * matmul (conjg(SMAT),DES)

! ... diagonals will be effectively real, but force it anyway:

      Do i=1,nopen; Q(i,i)=cmplx(real(Q(i,i)),0.0D0);  End do

! ... diagonalize the time-delay matrix:

      Call LAP_ZHEEV('V',nopen,nopen,Q,qval,info)

! ... time-delay is equal to largest eigenvalue:

      tdvalue = qval(nopen)

! ... decay mode is defined by corr. eigenvector:

      Do i = 1,nopen;  parg(i) = abs(Q(i,nopen))**2; End do

      Deallocate (SMAT1, SMAT2, SMAT, DES, Q, qval)

      End Subroutine time_del

