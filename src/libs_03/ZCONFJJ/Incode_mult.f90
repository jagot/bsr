!======================================================================
      Integer Function Incode_mult (itype,i1,i2)
!======================================================================
!     incode the integral
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: itype,i1,i2
      Integer, parameter :: jb = 2**10

      if(max(i1,i2).ge.jb) then
       write(*,'(a,2i5,a,i3)')  &
               ' Incode_mult: i1,i2 =',i1,i2,'  > base =',jb
       Stop
      end if
      if(itype.gt.3) then
       write(*,'(a,i5,a)')  &
               ' Incode_mult: itype =',itype,'  is out of range (3)'
       Stop
      end if

      Incode_mult = (itype*jb+i1)*jb+i2

      End Function Incode_mult


!======================================================================
      Subroutine Decode_mult (itype,i1,i2,int)
!======================================================================
!     decode the integral
!----------------------------------------------------------------------
      Implicit none
      Integer, parameter :: jb = 2**10
      Integer, intent(in) :: int
      Integer, intent(out) :: itype,i1,i2
      Integer :: k

      k  = int
      i2 = mod(k,jb);  k = k/jb
      i1 = mod(k,jb);  k = k/jb
      itype = mod(k,jb)

      End Subroutine Decode_mult
