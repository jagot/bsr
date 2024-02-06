!======================================================================
      Subroutine UPDATE_HS (ms,hm,ii,jj,ns,ks,d,sym)
!======================================================================
!     update 'small' channel block (ns,ns) in "big" channel block (ms,ms)
!     (ii,jj) = (1,1),(1,2),(2,1),(2,2) - define the "small" block
!     Possible array representation for input d-matrix:
!     sym = 's'  -->  symmetric banded upper-column storage mode
!     sym = 'n'  -->  non-symmetric band matrix
!     sym = 'x'  -->  non-symmetric full matrix
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: ii,jj,ms,ns,ks
      Real(8), intent(in) :: d(ns,*)
      Real(8), intent(inout) :: hm(ms,ms)
      Character(*), intent(in) :: sym
      Integer :: i,j, jp, imin,imax, i1,i2,j1,j2
      Real(8) :: y(ns,ns), yy(ns,ns)

      Select case(sym)

      Case('s')

       y = 0.d0
       do jp=1,ks;  do i=1,ns-jp+1;  j=i+jp-1
        y(i,j) = d(i,jp)
        y(j,i) = d(i,jp)
       end do; end do

      Case('n')

       y = 0.d0
       do jp = 1,ks+ks-1
         imin=max0( 1, 1 + ks-jp)
         imax=min0(ns,ns + ks-jp)
         do i = imin,imax;  j=i+jp-ks
           y(i,j) = d(i,jp)
         end do
       end do

      Case('x')

       y(1:ns,1:ns) = d(1:ns,1:ns)

      End Select

      Select Case (10*ii+jj)
       Case(11);  i1=1; i2=ns; j1=1; j2=ns
       Case(12);  i1=1; i2=ns; j1=ns+1;j2=ms
       Case(21);  i1=ns+1;i2=ms; j1=1; j2=ns
       Case(22);  i1=ns+1;i2=ms; j1=ns+1;j2=ms
       Case default; Stop 'unknown case in update_hs'
      End Select

      y = y/2.d0; yy = TRANSPOSE(y)
      hm(i1:i2,j1:j2) = hm(i1:i2,j1:j2) + y
      hm(j1:j2,i1:i2) = hm(j1:j2,i1:i2) + yy

      End Subroutine UPDATE_HS
