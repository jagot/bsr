!======================================================================
    Subroutine ZINTYm (nv,ks,ks1,ks2,B1,B2,ygr,mm,ym)
!======================================================================
!   Computes the array  ym(i,j) =  <B1_i | y(r) | B2 _j>
!
!   B1, B2 - two bases, spesified in gausian points ks
!            which may be different from ks1 or ks2
!
!   ygr - array of values of a specific function  y(r) at the gaussian
!         points of each interval, weighted by the gaussian weights
!----------------------------------------------------------------------
    Implicit none
    Integer, intent(in)  :: nv,ks,ks1,ks2,mm
    Real(8), intent(in)  :: B1(nv+1,ks,ks1),B2(nv+1,ks,ks2),ygr(nv,ks)
    Real(8), intent(out) :: ym(mm,mm)
    Integer :: i,j, iv, ith, jth

    ym = 0.d0
    Do iv = 1,nv                      ! over intervals
     Do ith = 1,ks1; i = iv+ith-1     ! over B1 splines
     Do jth = 1,ks2; j = iv+jth-1     ! over B2 splines
        ym(i,j) = ym(i,j) + SUM(ygr(iv,:)*b1(iv,:,ith)*b2(iv,:,jth))
     End do
     End do
    End do

    End Subroutine ZINTYM


!======================================================================
    Subroutine ZINTYk (k,ks1,ks2,B1,B2,mm,ym)
!======================================================================
!   Computes the array elements   <B1_i| r**k | B2_j>
!
!   B1, B2 - two bases, spesified in gausian points
!----------------------------------------------------------------------
    Use DBS_grid,  only: nv,ks
    Use DBS_gauss, only: gr, grw, gx, gw

    Implicit none
    Integer, intent(in)  :: k,ks1,ks2,mm
    Real(8), intent(in)  :: B1(nv+1,ks,ks1),B2(nv+1,ks,ks2)
    Real(8), intent(out) :: ym(mm,mm)

    Integer :: i,j, iv, ith,jth

    ym = 0.d0
    Do iv = 1,nv
     gx(:) = grw(iv,:) * gr(iv,:) ** k
     Do ith = 1,ks1; i = iv+ith-1
      gw (:) = gx(:) * b1(iv,:,ith)
     Do jth = 1,ks2; j = iv+jth-1
        ym(i,j) = ym(i,j) + SUM(gw(:)*b2(iv,:,jth))
     End do
     End do
    End Do

    End Subroutine ZINTYk




