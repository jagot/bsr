!=========================================================================
       Real(8) Function Coulomb_phase(l,z,k)
!=========================================================================
!      Coulomb phase is defined as:    arg Gamma[l+1+i(z/k)],  where
!       l - angular momentum
!       z - nuclear charge
!       k - wave number
!      We usually need the difference two Coulomb phases, so let define 
!      Coulomb phase shift as difference from sigma(0). Then
!      sigma(l+1)=sigma(l) + tg^-1( (z/k) / (l+1))  with sigma(0)=0
!-------------------------------------------------------------------------
       Implicit none
       Integer, intent(in) :: l
       Real(8), intent(in) :: z,k
       Integer :: i
       Real(8) :: q,s

       Coulomb_phase = 0.d0 
       if(l.le.0) Return
       if(z.eq.0.d0) Return         
       q = z/k
       Do i=1,l; s=DFLOAT(i)       
        Coulomb_phase = Coulomb_phase + DATAN2(q,s)
       End do
       
       End Function Coulomb_phase  !  do not checked yet

