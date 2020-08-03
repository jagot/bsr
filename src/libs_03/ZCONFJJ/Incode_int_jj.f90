!======================================================================
      Integer Function Incode_int(met,k,I1,I2,I3,I4)
!======================================================================
!     incode integral
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: met,k,I1,I2,I3,I4
      Integer :: ib2 = 2**2, ib5 = 2**5, ib10= 2**10

      if(max(i1,i2,i3,i4).ge.ib5) Stop 'Incode_int: i > pack-base'

      Incode_INT = ((((i1*ib5+i2)*ib5+i3)*ib5+i4)*ib10+k)*ib2+met

      End Function Incode_int


!======================================================================
      Subroutine Decode_int(met,k,I1,I2,I3,I4,int)
!======================================================================
!     decode the integral form "int_bnk" (relativistic)
!----------------------------------------------------------------------
      Implicit none
      Integer :: int, met,k,I1,I2,I3,I4, ii
      Integer :: ib2 = 2**2, ib5 = 2**5, ib10= 2**10

      ii = int
      met = mod(ii,ib2);  ii = ii/ib2
      k   = mod(ii,ib10); ii = ii/ib10
      I4  = mod(ii,ib5);  ii = ii/ib5
      I3  = mod(ii,ib5);  ii = ii/ib5
      I2  = mod(ii,ib5);  I1 = ii/ib5

      End Subroutine Decode_int

