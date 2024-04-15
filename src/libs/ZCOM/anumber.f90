!======================================================================
      Subroutine Anumber(number,na,A)
!======================================================================
!     Generate the character representation for digital number
!     (old version for portability, instead just use In.m format)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: number,na
      Character(na), intent(out) :: A
      Integer :: i,j,k,m,ii

      if(number.lt.0) Stop 'Anumber: number < 0 '

      ii = number
      Do i = na,1,-1
       k = 10**(i-1)
       j = ii/k
       if(j.gt.9) Stop 'Anumber: number too large'
       ii = ii - j*k
       m = na-i+1
       write(A(m:m),'(i1)') j
      End do

      End Subroutine Anumber
