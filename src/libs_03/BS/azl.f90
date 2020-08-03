!=======================================================================
    Real(8) Function azl(z,h,ks,lp)
!=======================================================================
!   Value of B_(lp+1)/r^lp at r = 0 where lp = l+1
!----------------------------------------------------------------------
    Implicit none
    Integer, intent(in) :: ks,lp
    Real(8), intent(in) :: z,h
    Integer :: j
    Real(8) :: c

    if (lp < ks ) then
      azl = 1.d0
      c = z/h
      Do j = 1,lp
        azl = (azl*c*(ks-j))/(j*j)
      End do
    else
      azl = 0.d0
    end if

    End Function azl
