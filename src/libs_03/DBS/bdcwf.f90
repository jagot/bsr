!=======================================================================
    Subroutine bdcwf_pq(n,k,z,Pcoef,Qcoef)
!=======================================================================
!   computes the spline expansion coefficients for hydrogenic radial
!   function for quantum numbers (n,k) and effective nuclear charge z.
!
!   Bases for P and Q may be different.
!
!   SUBROUTINE(S) called:  dcwf, gaussj
!-----------------------------------------------------------------------
!   on entry
!   --------
!       n,k     orbital quantum numbers
!       z       effective nuclear charge
!   on exit
!   -------
!       Pcoef   the spline coefficients for large component P(nl;r)
!       Qcoef   the spline coefficients for small component Q(nl;r)
!-----------------------------------------------------------------------
    Use DBS_grid,  only: ns,ks,nv, nsp,ksp, nsq,ksq
    Use DBS_gauss, only: gr,grw, pbsp,qbsp, fpbs,fqbs

    Implicit none
    Integer, intent(in)  :: n, k
    Real(8), intent(in)  :: z
    Real(8), intent(out) :: Pcoef(ns),Qcoef(ns)
    ! .. Local variables
    Integer :: i,iv,ip,m
    Real(8) :: E,pr(nv,ks),qr(nv,ks),a(ns,ns)

    ! .. obtain the values of the radial function
    ! .. at all the gaussian points
    ! .. and multiply by the gaussian weights

    Do m=1,ks
     Call dcwf(n,k,z,E,nv,gr(1,m),pr(1,m),qr(1,m))
     pr(:,m) = pr(:,m)*grw(:,m)
     qr(:,m) = qr(:,m)*grw(:,m)
    End do

    ! .. form the vector of inner products of the radial function and the
    ! .. spline basis functions

    Pcoef = 0.d0
    Do iv = 1,nv; Do ip = 1,ksp; i = iv+ip-1
     Pcoef(i) = Pcoef(i) + SUM(pr(iv,:)*pbsp(iv,:,ip))
    End do; End do

    Qcoef = 0.d0
    Do iv = 1,nv; Do ip = 1,ksq; i = iv+ip-1
     Qcoef(i) = Qcoef(i) + SUM(qr(iv,:)*qbsp(iv,:,ip))
    End do; End do

    ! .. solve the system of equations for coef

    a(1:nsp-1,1:nsp-1)=fpbs(2:nsp,2:nsp)
    Call gaussj (a,nsp-1,ns,Pcoef(2),1,1)  ! 1->ns  ?
    Pcoef(1)=0.d0

    a(1:nsq-1,1:nsq-1)=fqbs(2:nsq,2:nsq)
    Call gaussj (a,nsq-1,ns,Qcoef(2),1,1)  ! 1->ns  ?
    Qcoef(1)=0.d0

    End Subroutine bdcwf_pq

