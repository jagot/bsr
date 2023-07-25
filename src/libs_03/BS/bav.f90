!=======================================================================
      Subroutine bav(ns,ks,b,v,y,type,side)
!=======================================================================
!
!     Computes  y = B * v    or   y = v * B     (side =  'r' or 'l')
!
!     where b is a symmetric or non-symmetric banded matrix
!     (type = 's' or 'n')
!
!     and  v, y  are vectors
!
!-----------------------------------------------------------------------
!
!     on entry
!     --------
!       ks      the number of diagonals
!       ns      the order of the matrix
!       B       the  banded matrix in column (upper) storage mode
!       v       vector
!
!     on exit
!     -------
!       y       y = B*v
!
!-----------------------------------------------------------------------

      Implicit none
      Integer, intent(in) :: ns, ks
      REAL(8), intent(in) :: v(ns)
      REAL(8), intent(in) :: b(ns,*)
      REAL(8), intent(out) :: y(ns)
      Character :: type, side

      ! .. Local variables

      Integer :: i, j, ip, imin, imax


      if(type.eq.'s') then

        do i=1,ns
          y(i) = b(i,1)*v(i)
        end do

        do ip = 2,ks
          do i = 1,ns-ip+1
            j = i+ip-1
            y(i) = y(i) + b(i,ip)*v(j)
            y(j) = y(j) + b(i,ip)*v(i)
          end do
        end do

      elseif(type.eq.'n'.and.side.eq.'r') then

        y = 0.d0
        do ip = 1,ks+ks-1
          imin=max( 1, 1 + ks-ip)
          imax=min(ns,ns + ks-ip)
          do i = imin,imax
            j = i+ip-ks
            y(i) = y(i) + b(i,ip)*v(j)
          end do
        end do

      elseif(type.eq.'n'.and.side.eq.'l') then

        y = 0.d0
        do ip = 1,ks+ks-1
          imin=max( 1, 1 + ks-ip)
          imax=min(ns,ns + ks-ip)
          do i = imin,imax
            j = i+ip-ks
            y(j) = y(j) + b(i,ip)*v(i)
          end do
        end do

      end if

      End Subroutine bav
