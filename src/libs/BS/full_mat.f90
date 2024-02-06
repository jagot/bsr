!======================================================================
      Subroutine Full_mat_sym(ns,ks,x,y,sym)
!======================================================================
! ... provide full matrix from band format
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ns,ks
      REAL(8) :: x(ns,*),y(ns,ns)
      CHARACTER(1) :: sym
      INTEGER ::  i,j, jp, imin,imax

      y = 0.d0

      Select case(sym)

      Case('s')       ! symetric, column

       do jp=1,ks;  do i=1,ns-jp+1;  j=i+jp-1
        y(i,j) = x(i,jp)
        y(j,i) = x(i,jp)
       end do; end do

      Case('l')       ! symetric, row

       do jp=1,ks;  do i=ks+1-jp,ns;  j=i+jp-ks
        y(i,j) = x(i,jp)
        y(j,i) = x(i,jp)
       end do; end do

      Case('n')       ! nonsymetric

       do jp = 1,ks+ks-1
         imin=max0( 1, 1 + ks-jp)
         imax=min0(ns,ns + ks-jp)
         do i = imin,imax;  j=i+jp-ks
           y(i,j) = x(i,jp)
         end do
       end do

      Case('x')       !  full

        y(1:ns,1:ns) = x(1:ns,1:ns)

      End Select

      End Subroutine Full_mat_sym
