!====================================================================
      Subroutine NUM (z,k1,k2,i_max,acr)
!====================================================================
!    finds integer numbers k1 and k2 such that abs(z) = sqrt(k1/k2)
!    and sign(z) = sign(k1)
!    k1, k2 < i_max  - restrictions on possible values of k1,k2
!    acr - tollerance for above fitting
!--------------------------------------------------------------------
      Implicit none
      Real(8), intent(in)  :: z,acr
      Integer, intent(in)  :: i_max
      Integer, intent(out) :: k1,k2
      Integer :: i_min, i1,i2, m1,m2
      Real(8) :: zz, s,ss, s1,s2

      k1=0;  k2=1; zz=z*z
      if(zz.lt.1.0/i_max) Return

      s=zz
      i_min=zz
      if(i_min.lt.1) i_min=1
      Do i1=i_min,i_max
       s1=DBLE(i1)
       m1=s1/zz;  if(m1.lt.1) m1=1;  m2=m1+1
      Do i2=m1,m2
       s2=DBLE(i2)
       ss=abs(s1/s2-zz)
       if(ss.lt.s) then
        s = ss
        k1 = i1; if(z.lt.0.d0) k1 = -k1;   k2 = i2
       end if
       if(ss.lt.acr) Return
      End do
      End do

      End Subroutine NUM
