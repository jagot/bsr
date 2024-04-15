!======================================================================
      Subroutine Sym_mat (n,m,R,dmn,dmx,dav)
!======================================================================
!     evaluate symmetry of the matrix R(n,n), and then symmetrize it
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,m
      Real(8), intent(inout) :: R(m,m)
      Real(8), intent(out) :: dmn,dmx,dav
      Integer :: i,j,k
      Real(8) :: dij

      dmx=0.d0; dmn=0.d0; dav=0.d0; if(n.le.1) Return

      k = 0
      Do i = 2,n
       Do j = 1,i-1
        dij = R(i,j)+R(j,i); if(dij.eq.0.d0) Cycle
        dij = dabs ((R(i,j)-R(j,i))/dij)
        dav = dav + dij; k=k+1
        if(dmx.eq.0.d0.or.dij.gt.dmx) dmx = dij
        if(dmn.eq.0.d0.or.dij.lt.dmn) dmn = dij
       End do
      End do

      if(k.gt.0) dav = dav / k
      if(dmx.gt.1.d-15) dmx = LOG10(dmx/2.d0)
      if(dmn.gt.1.d-15) dmn = LOG10(dmn/2.d0)
      if(dav.gt.1.d-15) dav = LOG10(dav/2.d0)

      Do i = 1,n-1
       Do j = i+1,n
        R(i,j)=(R(i,j)+R(j,i))/2.d0;  R(j,i) = R(i,j)
       End do
      End do

      End Subroutine Sym_mat