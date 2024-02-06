!======================================================================
      Subroutine Clean_a(a)
!======================================================================
!     remove blanks from the character string
!----------------------------------------------------------------------
      Character(*) :: a

      len = len_trim(a)
      m = 0
      Do i=1,len
       if(a(i:i).eq.' ') cycle
       m = m + 1
       a(m:m)=a(i:i)
      End do
      if(m.lt.len) a(m+1:len) = ' '

      End Subroutine Clean_a
