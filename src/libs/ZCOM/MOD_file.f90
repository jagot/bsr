!====================================================================
      Module internal_file
!====================================================================
!     contains a set of character lines
!--------------------------------------------------------------------
      Implicit none 
    
      Integer :: mlines =  0      !  max. number of lines
      Integer :: nlines =  0      !  current number of lines
      Integer :: klines = 2**8    !  initial space  
      Integer, parameter :: mlen = 2**7   !  max. length

      Character(mlen), allocatable :: aline(:)

      End Module internal_file


!======================================================================
      Subroutine alloc_file(m)
!======================================================================
!     allocate, deallocate or reallocate arrays in module "internal_file"
!----------------------------------------------------------------------
      Use internal_file

      Implicit none
      Integer, intent(in) :: m
      Character(mlen), allocatable :: bline(:)

      if(m.le.0) then
       if(allocated(aline)) Deallocate (aline)
       mlines = 0; nlines = 0; Return
      elseif(.not.allocated(aline)) then
       mlines = m; Allocate(aline(m))
      elseif(m.le.mlines) then
       Return
      elseif(nlines.eq.0) then
       Deallocate (aline);  mlines = m
       Allocate(aline(m))
      else
       Allocate(bline(nlines))
       bline(1:nlines)=aline(1:nlines); Deallocate(aline)
       Allocate(aline(m)); aline(1:nlines)=bline(1:nlines)
       Deallocate(bline)
       mlines = m
      end if

      End Subroutine alloc_file 


!=======================================================================
      Integer Function Iadd_line(line)
!=======================================================================
!     add the integral to the list 'internal_file'
!-----------------------------------------------------------------------
      Use internal_file

      Implicit none
      Character(*) :: line
      Integer :: i

      if(mlines.eq.0) Call Alloc_file(klines)

! ... check if the same integral is already in list:

      Do i=1,nlines
       if(aline(i).ne.line) Cycle
       Iadd_line=i; Return
      End do

! ... add new integral:

      if(nlines.eq.mlines) Call Alloc_file(mlines+klines)
      nlines=nlines+1; aline(nlines)=line
      Iadd_line = nlines

      End Function Iadd_line 
