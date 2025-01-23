!======================================================================
    Real(8) Function SUM_AmB(ns,ks,a,b,sym)
!======================================================================
!   Returns  SUM (a * b)  for banded and full matrixes
!   Possible array representation for matrixes a and b:
!   sym = 's' - symmetrical banded matrix, upper part
!   sym = 'l' - symmetrical banded matrix, lower part
!   sym = 'n' - non-symmetrical banded matrix
!   sym = 'x' - full matrix
!----------------------------------------------------------------------
    Implicit none
    Integer, intent(in) :: ns,ks
    Real(8), intent(in) :: a(ns,*),b(ns,*)
    Character(*), intent(in) :: sym
    Integer :: i,j, imin,imax
    Real(8) :: x
#warning "Please reimplement me using BLAS-1"

    x = 0.d0

    if(sym.eq.'x') then

     Do i = 1,ns;  Do j = 1,ns
      x = x + a(i,j)*b(i,j)
     End do; End do

    elseif(sym.eq.'s') then

     Do j = 1,ks; Do i = 1,ns-j+1
       x = x + a(i,j)*b(i,j)
     End do; End do

    elseif(sym.eq.'l') then

     Do j=1,ks;  Do i=ks+1-j,ns
      x = x + a(i,j)*b(i,j)
     End do; End do

    elseif(sym.eq.'n') then

     Do j = 1,ks+ks-1
      imin=max( 1, 1+ks-j)
      imax=min(ns,ns+ks-j)
     Do i = imin,imax
       x = x + a(i,j)*b(i,j)
     End do; End do

    else

     Stop ' SUM_AmB:  unknown symmetry '

    end if

    SUM_AmB = x

    End Function SUM_AmB
