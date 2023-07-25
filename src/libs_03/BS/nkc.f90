!======================================================================
      REAL(8) FUNCTION nkc (i1,j1,i2,j2,k)
!======================================================================
!               k
!     Returns  N (i1, j1; i2, j2) - direct summations of moments
!                                   over cells 
!----------------------------------------------------------------------

      USE spline_param
      USE spline_orbitals, p => pbs
      USE spline_moments
      USE spline_atomic
    
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i1,j1,i2,j2,k
    
      ! .. local variables
    
      INTEGER(4) :: i,j, iv, ii, ivi,ivj, ik
      REAL(8), DIMENSION(nv) :: v1, v2
      REAL(8), DIMENSION(ks*(ks+1)/2) :: a,b
      REAL(8) :: s
    
      ! .. moments calculations
    
      Call nk_moments(k)
    
      ik=ks*(ks+1)/2
      
      nkc = 0.d0
    
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
    
       v1(iv) = SUM(a*rkd2(1:ik,iv))
       v2(iv) = SUM(b*rkd1(1:ik,iv))
    
       ! .. diagonal contributions
    
       Do j=1,ik
        nkc = nkc + SUM(a(1:ik)*rkd(1:ik,j,iv))*b(j)
       End do
    
      End do
    
      ! .. the upper and lower regions
    
      s = 0.d0
      Do iv = 2,nv
       s = s + v2(iv-1)
       nkc = nkc + s*v1(iv)
      End do
    
      nkc = nkc * fine
    
      End FUNCTION nkc

