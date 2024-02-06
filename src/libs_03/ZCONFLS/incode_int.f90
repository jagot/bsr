!======================================================================
      Integer Function Incode_INT (met,k,I1,I2,I3,I4)
!======================================================================
!     incode the two-electron integral:  Rk(i1,i2;i3,i4)
!     met -> type of integral
!----------------------------------------------------------------------
      Implicit none
      Integer :: ib4=2**4, ib5=2**5, ib6=2**6
      Integer, intent(in) :: met,k,I1,I2,I3,I4
      if(max(i1,i2,i3,i4).ge.ib5) Stop 'Incode_INT: i >= 2**5'
      if(met.gt.ib4) Stop 'Incode_INT: met >= 2**4'
      if(k.gt.ib6) Stop 'Incode_INT: met >= 2**6'

      Incode_INT = ((((k*ib4+met)*ib5+i4)*ib5+i3)*ib5+i2)*ib5+i1

      End Function Incode_INT

!======================================================================
      Subroutine Decode_INT (met,k,I1,I2,I3,I4,int)
!======================================================================
!     decode the two-electron integral:  Rk(i1,i2;i3,i4)
!     met -> type of integral
!----------------------------------------------------------------------
      Implicit none
      Integer :: ib4=2**4, ib5=2**5
      Integer, Intent(in)  :: int
      Integer, Intent(out) :: met,k,I1,I2,I3,I4

      k   = int
      I1  = mod(k,ib5);  k = k/ib5
      I2  = mod(k,ib5);  k = k/ib5
      I3  = mod(k,ib5);  k = k/ib5
      I4  = mod(k,ib5);  k = k/ib5
      met = mod(k,ib4);  k = k/ib4

      End Subroutine Decode_INT

!======================================================================
      Subroutine Decode_met(met,int)
!======================================================================
!     decode just "type" of the two-electron integral
!----------------------------------------------------------------------
      Implicit none
      Integer :: ib4=2**4, ib20=2**20
      Integer, Intent(in)  :: int
      Integer, Intent(out) :: met

      met = int/ib20; met = mod(met,ib4)

      End Subroutine Decode_met

!======================================================================
      Integer Function Incode_mult (itype,i1,i2)
!======================================================================
!     incode the multipole integral:  d (i1,i2)
!     itype = 0  -  overlap
!           = 1  -  electric multipole
!           = 2  -  magnetic multipole
!           = 3  -  ???
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: itype,i1,i2
      Integer :: jb=2**10

      if(max(i1,i2).ge.jb) then
       write(*,'(a,2i5,a,i3)')  &
               ' Incode_mult: i1,i2 =',i1,i2,'  > base =',jb
       Stop
      end if
      if(itype.gt.3) then
       write(*,'(a,i5,a)')  &
               ' Incode_mult: itype =',itype,'  is out of range (>3)'
       Stop
      end if

      Incode_mult = (itype*jb+i1)*jb+i2

      End Function Incode_mult


!======================================================================
      Subroutine Decode_mult (itype,i1,i2,int)
!======================================================================
!     decode the multipole integral:  d (i1,i2)
!     itype = 0  -  overlap
!           = 1  -  electric multipole
!           = 2  -  magnetic multipole
!           = 3  -  ???
!----------------------------------------------------------------------
      Implicit none
      Integer, Intent(in) :: int
      Integer, Intent(out) :: itype,i1,i2
      Integer :: jb=2**10, k

      k  = int
      i2 = mod(k,jb);  k = k/jb
      i1 = mod(k,jb);  k = k/jb
      itype = mod(k,jb)

      End Subroutine Decode_mult



