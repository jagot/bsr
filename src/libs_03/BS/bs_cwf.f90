!======================================================================
      Subroutine bs_cwf(z,l,k,c,acc)
!======================================================================
!     B-spline expansion for regular Coulomb function with E = k^2
!----------------------------------------------------------------------
!     z - charge
!     l - orbital momentum
!     k - linear momemtum
!     c - output expansion vector
!     acc - accuracy
!     Normalization and estimation of accuracy is given by comparison
!     with COULFG at last interval
!----------------------------------------------------------------------
      Use spline_param
      Use spline_grid
      Use spline_galerkin, bbb => bb

      Implicit none
      Integer, intent(in) :: l
      Real(8), intent(in) :: z,k
      Real(8), intent(out) :: acc,c(ns)

      Integer :: i,j,jp,m,n,ith,ifail,info,ii
      Integer, allocatable  :: IPIV(:)
      Real(8) :: fl,eta,rho,cn,acc1,acc2,ee,zz,fc,gc,fcp,gcp
      REAL(8), allocatable :: a(:,:),b(:,:),aa(:,:),bb(:,:),hl(:,:)
      REAL(8), allocatable :: cf(:),cfb(:),cc(:),c1(:),c2(:)

      Allocate(a(ns,ns),b(ns,ns),aa(ns,ns),bb(ns,ns),hl(ns,ks))
      Allocate(IPIV(ns),cf(ks),cfb(ks),cc(ns),c1(ns),c2(ns))

! ... find reper points:

      fl = DBLE(l)
      eta = -z/k
      Do i = 1,ks
       rho = k*gr(nv,i)
       Call Zcoulfg90 (rho,eta,fl,fc,gc,fcp,gcp,0,ifail)
       cf(i) = fc
      End do

! ... set up hl matrix

      ee = k * k /2.d0
      fl = l*(l+1.d0)
      zz = 2.d0*z
      hl = db2 - fl*rm2 + zz*rm1

! ... store the (ns-1,ns) asymmetric value in hl(1,1)

      hl(1,1) = hl(ns,ks-1) + (db2(1,1)-db2(ns,ks-1))

! ... atomic units:

      hl = -0.5d0 * hl

! ... full matrixes:

      a = 0.d0
      b = 0.d0
      Do j = 1,ks;  Do i = ks-j+1,ns; jp=i-ks+j
        a(i,jp)=hl(i,j);  a(jp,i)=hl(i,j)
        b(i,jp)=sb(i,j);  b(jp,i)=sb(i,j)
      End do; End do
      a(ns-1,ns)=hl(1,1)

! ... delete first B-splines:

      ii = l+1; if(ii.ge.ks) ii=1
      n = ns - ii

! ... first option:  c(n-1) = 1

      aa(1:n,1:n) = a(ii+1:ns,ii+1:ns)
      bb(1:n,1:n) = b(ii+1:ns,ii+1:ns)
      aa(n-1,n) = aa(n,n-1)
      aa = aa - ee*bb
      cc = 0.d0;  cc(n-1) = 1.d0
      Call DGESV( n,1, aa, ns, IPIV, cc, ns, INFO )
      c1 = 0.d0; c1(ii+1:ns) = cc(1:n)
      acc1 = 1.d0
      if(INFO.eq.0) then
      cfb = 0.d0
      Do m = 1,ks                 ! over gausian points
       Do ith = 1,ks              ! over B-splines in given interval
        i = nv+ith-1              ! B-spline index
        cfb(m) = cfb(m) + c1(i)*bsp(nv,m,ith)
        End do
      End do
      cn = 0.d0
      Do m = 1,ks; cfb(m) = cf(m)/cfb(m); cn=cn+cfb(m); End do
      cn=cn/ks; c1 = cn*c1
      acc1=0.d0
      Do m=1,ks; acc1=acc1 + abs(cn-cfb(m)); End do
      acc1=acc1/cn/ks; acc1=abs(acc1)
      end if

! ... second option:  c(n) = 1

      aa(1:n,1:n) = a(ii+1:ns,ii+1:ns)
      bb(1:n,1:n) = b(ii+1:ns,ii+1:ns)
      aa(n,n-1) = aa(n-1,n)
      aa = aa - ee*bb
      cc = 0.d0; cc(n) = 1.d0
      Call DGESV( n,1, aa, ns, IPIV, cc, ns, INFO )
      c2 = 0.d0; c2(ii+1:ns) = cc(1:n)
      acc2 = 1.d0
      if(INFO.eq.0) then
      cfb = 0.d0
      Do m = 1,ks                 ! over gausian points
       Do ith = 1,ks              ! over B-splines in given interval
        i = nv+ith-1              ! B-spline index
        cfb(m) = cfb(m) + c2(i)*bsp(nv,m,ith)
        End do
      End do
      cn = 0.d0
      Do m = 1,ks; cfb(m) = cf(m)/cfb(m); cn=cn+cfb(m); End do
      cn=cn/ks; c2 = cn*c2
      acc2=0.d0
      Do m=1,ks; acc2=acc2 + abs(cn-cfb(m)); End do
      acc2=acc2/cn/ks; acc2=abs(acc2)
      end if

      if(acc1.gt.acc2) then
       acc = acc2; c = c2
      else
       acc = acc1; c = c1
      end if

      End Subroutine bs_cwf

