!     don't confuse bav and bxv - different storage modes  !!!

!=======================================================================
      Subroutine bav(ns,ks,B,v,y,type,side)
!=======================================================================
!     Computes  y = B * v    or   y = v * B     (side =  'r' or 'l')
!
!     where B is a symmetric or non-symmetric banded matrix
!     (type = 's' or 'n') in upper-band storage mode
!     and  v, y  are vectors      (used for RK-integrals)
!-----------------------------------------------------------------------
      Implicit none

      Integer, intent(in)  :: ns,ks
      Real(8), intent(in)  :: v(ns)
      Real(8), intent(out) :: y(ns)
      Real(8), intent(in)  :: b(ns,*)
      Character(1) :: type, side

      Integer :: i, j, ip, imin, imax

      if(type.eq.'s') then

       y(1:ns) = b(1:ns,1)*v(1:ns)

       do ip = 2,ks
        do i = 1,ns-ip+1;  j = i+ip-1
          y(i) = y(i) + b(i,ip)*v(j)
          y(j) = y(j) + b(i,ip)*v(i)
        end do
       end do

      elseif(side.eq.'r') then

        y = 0.d0
        do ip = 1,ks+ks-1
          imin=max( 1, 1 + ks-ip)
          imax=min(ns,ns + ks-ip)
          do i = imin,imax;  j = i+ip-ks
            y(i) = y(i) + b(i,ip)*v(j)
          end do
        end do

      elseif(side.eq.'l') then

        y = 0.d0
        do ip = 1,ks+ks-1
          imin=max( 1, 1 + ks-ip)
          imax=min(ns,ns + ks-ip)
          do i = imin,imax; j = i+ip-ks
            y(j) = y(j) + b(i,ip)*v(i)
          end do
        end do

      end if

      End Subroutine bav


!=======================================================================
      Subroutine bxv(k,n,b,v,y)
!=======================================================================
!     Computes   y = b * v    where b is a symmetric, banded matrix,
!     in lower-band storage mode,  and v, y are vectors.
!     (used for L-integrals)
!-----------------------------------------------------------------------
      Implicit none 
      Integer, intent(in)  :: n,k
      Real(8), intent(in)  :: v(n),b(n,k)
      Real(8), intent(out) :: y(n)
      Integer :: i, j, jp

! ... contribution from central diagonal (jp=k)

      Do i=1,n;  y(i) = b(i,k)*v(i) ; End do

! ... off diagonal

      Do jp = 1,k-1
       Do i = k-jp+1,n
        j = i-k+jp
        y(i) = y(i) + b(i,jp)*v(j)             ! sub_diagonals
        y(j) = y(j) + b(i,jp)*v(i)             ! super-diagonals
       End do
      End do

      End Subroutine bxv

