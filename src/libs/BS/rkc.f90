!======================================================================
      REAL(8) FUNCTION rkc (i1,j1,i2,j2,k)
!======================================================================
!               k
!     Returns  R (i1, j1; i2, j2) - direct summation of moments
!                                   over cells
!----------------------------------------------------------------------

      USE spline_param
      USE spline_orbitals, p => pbs
      USE spline_moments

      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: i1,j1,i2,j2,k

      ! .. local variables

      INTEGER(4) :: i,j, iv, ii, ivi,ivj, ik
      REAL(8), DIMENSION(nv) :: v1,v2,v3,v4
      REAL(8), DIMENSION(ks*(ks+1)/2) :: a,b
      REAL(8) :: s1, s2

      !    moments calculations

      Call rk_moments(k)

      ik = ks*(ks+1)/2

      rkc= 0.d0

      Do iv = 1,nv

       ii=0
       Do i=1,ks
        ivi=iv+i-1
        Do j=i,ks
         ivj=iv+j-1
         ii = ii+1
         if(i.eq.j) then
          a(ii) = p(ivi,i1)*p(ivj,i2)
          b(ii) = p(ivi,j1)*p(ivj,j2)
         else
          a(ii) = p(ivi,i1)*p(ivj,i2) + p(ivi,i2)*p(ivj,i1)
          b(ii) = p(ivi,j1)*p(ivj,j2) + p(ivi,j2)*p(ivj,j1)
         end if
        End do
       End do

       v1(iv) = SUM(a*rkd1(1:ik,iv))
       v2(iv) = SUM(b*rkd1(1:ik,iv))
       v3(iv) = SUM(a*rkd2(1:ik,iv))
       v4(iv) = SUM(b*rkd2(1:ik,iv))

       ! .. diagonal contributions

       Do j=1,ik
        rkc = rkc + SUM(a(1:ik)*rkd(1:ik,j,iv))*b(j)
       End do

      End do

      ! the upper and lower regions

      s1 = 0.d0
      s2 = 0.d0
      Do iv =  2,nv
        s1 = s1 + v1(iv-1)
        rkc = rkc + s1*v4(iv)
        s2 = s2 + v2(iv-1)
        rkc = rkc + s2*v3(iv)
      End do

      End FUNCTION rkc

