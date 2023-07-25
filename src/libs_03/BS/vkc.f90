!======================================================================
      Real(8) Function vkc(i1,j1,i2,j2,k)
!======================================================================
!               k
!     Returns  V (i1, j1; i2, j2) - direct summations of moments
!                                   over cells 
!----------------------------------------------------------------------
      USE spline_param
      USE spline_orbitals, p => pbs
      USE spline_moments
      USE spline_atomic
   
      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
   
      ! .. local variables
   
      INTEGER :: i,j, iv, ivi, ivj, ii, ik,jk 
      REAL(8), DIMENSION(nv) :: v1, v2, v3, v4
      REAL(8), DIMENSION(ks*ks) ::  a, b
      REAL(8) :: s1, s2
   
      ! .. check the need of calculations
                  
      Call vk_moments(k)
   
      ik = ks*(ks+1)/2
      jk = ks*ks
   
      vkc = 0.d0
   
      Do iv = 1,nv
      
       ii=0
       Do i=1,ks
        ivi=iv+i-1
        Do j=1,ks
         ivj=iv+j-1
         ii = ii+1
         a(ii) = p(ivi,i1)*p(ivj,i2)
        End do
       End do
   
       ii=0
       Do i=1,ks
        ivi=iv+i-1
        Do j=i,ks
         ivj=iv+j-1
         ii = ii+1
         if(i.eq.j) then
          b(ii) = p(ivi,j1)*p(ivj,j2)
         else
          b(ii) = p(ivi,j1)*p(ivj,j2) + p(ivi,j2)*p(ivj,j1)
         end if
        End do
       End do
   
       v1(iv) = SUM(rkd3(1:jk,iv)*a)
       v2(iv) = SUM(rkd2(1:ik,iv)*b)
       v3(iv) = SUM(rkd4(1:jk,iv)*a)
       v4(iv) = SUM(rkd1(1:ik,iv)*b)
       
       ! the diagonal cell contribution
         
       Do j = 1,ik
        vkc = vkc + SUM(a(1:jk)*rkd(1:jk,j,iv))*b(j)
       End do
       
      End do 
       
      ! the upper and lower regions
   
      s1 = 0.d0
      s2 = 0.d0
      Do iv =  2,nv
       s1 = s1 + v1(iv-1)
       vkc = vkc + s1*v2(iv)
       s2 = s2 + v4(iv-1)
       vkc = vkc + s2*v3(iv)
      End do     
       
      vkc = vkc * fine
   
      End Function vkc

