!=======================================================================
      Real(8) FUNCTION grad(i,j)
!=======================================================================!
!     <p(i)| d + [l(j)(l(j)+1)-l(i)*(l(i)+1)]/2r |p(j)>
!-----------------------------------------------------------------------
      USE spline_param
      USE spline_galerkin
      USE spline_orbitals, L => LBS, mx => MBS, p => PBS

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i,j
      INTEGER :: ll,ii,jj,jp

      if(iabs(L(I)-L(J)).ne.1)  Stop ' BS_GRAD: L(i) - L(j) <> 1'

      ll = (L(j)*(L(j)+1)-L(i)*(L(i)+1))/2

      grad = ll * SUM(p(:,i)*p(:,j)*rm1(:,ks))

      do jp=1,ks-1;  do ii=ks+1-jp,ns;  jj=jp+ii-ks
         grad =  grad &
           +     db1(ii,jp) * (p(jj,j)*p(ii,i) - p(jj,i)*p(ii,j)) &
           +  ll*rm1(ii,jp) * (p(jj,j)*p(ii,i) + p(jj,i)*p(ii,j))
      end do;  end do

      END FUNCTION grad
