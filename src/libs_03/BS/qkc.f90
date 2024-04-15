!======================================================================
      Real(8) FUNCTION qkc(i1,j1,i2,j2,k)
!======================================================================
!                 k
!     Return  Q (i1, j1; i2, j2) - direct summation of moments
!                                  over shell
!
!     Calls:  qk_moments          non-correct yet!!!
!----------------------------------------------------------------------

      USE spline_param
      USE spline_orbitals, p => pbs
      USE spline_moments
      USE spline_atomic

      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: i1,j1,i2,j2,k

      ! .. local variables

      INTEGER(4) :: i,j, iv, ii, ivi,ivj, ik
      REAL(8), DIMENSION(nv) :: v1, v2, v3, v4
      REAL(8), DIMENSION(ks*ks) :: a,b
      REAL(8) :: s1, s2

      ! .. moments calculations

      Call qk_moments(k)

      ik=ks*ks

      qkc = 0.d0

      Do iv = 1,nv

       ii=0
       Do i=1,ks
        ivi=iv+i-1
        Do j=1,ks
         ivj=iv+j-1
         ii = ii+1
         a(ii) = p(ivi,i1)*p(ivj,i2)
         b(ii) = p(ivi,j1)*p(ivj,j2)
        End do
       End do

       v1 (iv) = SUM(a*rkd3(:,iv))
       v2 (iv) = SUM(b*rkd2(:,iv))
       v3 (iv) = SUM(a*rkd4(:,iv))
       v4 (iv) = SUM(b*rkd1(:,iv))

       ! .. diagonal cell contribution

       Do j=1,ik
         qkc = qkc + SUM(a(1:ik)*rkd(1:ik,j,iv))*b(j)
       End do

      End do

      ! the upper and lower regions

      s1 = 0.d0
      s2 = 0.d0
      Do iv =  2,nv
       s1 = s1 + v1(iv-1)
       qkc = qkc + s1*v2(iv)
       s2 = s2 + v4(iv-1)
       qkc = qkc + s2*v3(iv)
      End do

      qkc = qkc * fine

      End FUNCTION qkc

