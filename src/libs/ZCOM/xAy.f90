!=====================================================================
      Real(8) Function xAy (n,k,array,type,x,y)
!=====================================================================
!
!     Returns   <x | array |y>  for different array types
!
!     The original size of array is n*n. The bandwidth of array is 2k-1
!     for the non-symmetric case, and is k for symmetric case.
!     The lower-band column storage mode is assumed.
!---------------------------------------------------------------------
      Integer, intent(in) :: n, k
      Character, intent(in) :: type
      Real(8), intent(in) :: array(n,*), x(n), y(n)
      Integer :: i

      xAy = 0.d0

      if ( type .eq. 'f') then
       Do i=1,n
        xAy = xAy + SUM(array(1:n,i)*x(1:n)) * y(i)
       End do
      else
       Stop  'xAy: not yet supported array type'
      end if

      End Function xAy


!=====================================================================
      Real(8) Function  VAV (ni,ki,nj,kj,n,a,x,y,sym)
!=====================================================================
!            vav = <x |array | y>
!---------------------------------------------------------------------
!     n      the leading dimension for array "a"
!     ni,nj  number of splines
!     ki,kj  band width
!     a      banded matrix (in diff.storage modes):
!     sym    'x' general
!            'l' symmetric, lower-band storage mode
!            'u' symmetric, upper-band storage mode
!            'n' non-symmetric band storage mode
!     x      vector with length ni
!     y      vector with length nj
!----------------------------------------------------------------------

      Integer, intent(in)   :: n,ni,ki,nj,kj
      Real(8), intent(in)   :: a(n,*), x(*), y(*)
      Character, intent(in) :: sym

      Integer :: i,j,jp
      Real(8) :: v(n)

      if (sym.eq.'x') then

       Do j=1,nj
        v(j)=SUM(x(1:ni)*a(1:ni,j))
       End do
       vav = SUM(v(1:nj)*y(1:nj))

      elseif(sym.eq.'u') then

       vav = SUM(a(1:ni,1)*x(1:ni)*y(1:ni))
       Do jp=2,ki;  Do i=1,ni-jp+1; j=i+jp-1
        vav = vav + a(i,jp) * (x(i)*y(j) + x(j)*y(i))
       End do; End do

      elseif(sym.eq.'l') then

       vav = SUM(a(1:ni,ki)*x(1:ni)*y(1:ni))
       Do jp=1,ki-1;  Do i=ki+1-jp,ni; j=i+jp-ki
        vav = vav + a(i,jp) * (x(i)*y(j) + x(j)*y(i))
       End do; End do

      elseif(sym.eq.'n') then

       vav = 0.d0
       Do jp = 1,k+k-1
        i1 = max0(1,1+k-jp)
        i2 = min0(n,n+k-jp)
        Do i = i1,i2;  j = i+jp-k
         vav = vav + array(i,jp)*x(i)*y(j)
        End do
       End do

      else

        Stop 'VAV: unknown symmetry'

      end if

      End Function vav
