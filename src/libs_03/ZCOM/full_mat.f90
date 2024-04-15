!======================================================================
      Subroutine Full_mat_sym(ns,ks,x,y,sym)
!======================================================================
!     Builds full matrix from different storage modes for banded matrix
!     (see p.328 in BSR description, CPC 174 (2006))
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: ns,ks
      Real(8), intent(in) :: x(ns,*)
      Real(8), intent(out) :: y(ns,ns)
      Character(*), intent(in) :: sym
      Integer ::  i,j, jp, imin,imax

      y = 0.d0

      Select case(sym)

      Case('s')   ! symmetric upper-column storage

       do jp=1,ks;  do i=1,ns-jp+1;  j=i+jp-1
        y(i,j) = x(i,jp)
        y(j,i) = x(i,jp)
       end do; end do

      Case('l')   ! symmetric lower-column storage

       do jp=1,ks;  do i=ks+1-jp,ns;  j=i+jp-ks
        y(i,j) = x(i,jp)
        y(j,i) = x(i,jp)
       end do; end do

      Case('n')   !  non-symmetric banded

       do jp = 1,ks+ks-1
         imin=max0( 1, 1 + ks-jp)
         imax=min0(ns,ns + ks-jp)
         do i = imin,imax;  j=i+jp-ks
           y(i,j) = x(i,jp)
         end do
       end do

      Case('x')   ! just full (for completeness)

        y(1:ns,1:ns) = x(1:ns,1:ns)

      End Select

      End Subroutine Full_mat_sym
