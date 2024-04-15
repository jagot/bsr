!======================================================================
    REAL(8) FUNCTION SUM_AmB(ns,ks,a,b,sym)
!======================================================================
!   Returns  SUM ( a * b)  for banded and full matrixes
!
!   sym = 's' - symmetrical banded matrix (column upper-band storage)
!   sym = 'l' - symmetrical banded matrix (column lower-band storage)
!   sym = 'n' - non-symmetrical banded matrix
!   sym = 'x' - full matrix
!----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, INTENT(in) :: ns,ks
    CHARACTER(1), INTENT(in) :: sym
    REAL(8), INTENT(in), DIMENSION(ns,*) :: a,b

    INTEGER :: i,j, imin,imax
    REAL(8) :: x

    x = 0.d0

    if(sym.eq.'x') then

     Do i = 1,ns
      Do j = 1,ns
       x = x + a(i,j)*b(i,j)
      End do
     End do

    elseif(sym.eq.'s') then

     Do j = 1,ks
      Do i = 1,ns-j+1
       x = x + a(i,j)*b(i,j)
      End do
     End do

    elseif(sym.eq.'l') then

     Do j = 1,ks
      Do i = ks+1-j,ns
       x = x + a(i,j)*b(i,j)
      End do
     End do

    elseif(sym.eq.'n') then

     Do j = 1,ks+ks-1
      imin=max( 1, 1+ks-j)
      imax=min(ns,ns+ks-j)
      Do i = imin,imax
       x = x + a(i,j)*b(i,j)
      End do
     End do

    else

     Stop ' SUM_AmB:  unknown symmetry '

    end if

    SUM_AmB = x

    END FUNCTION SUM_AmB
