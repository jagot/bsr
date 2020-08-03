
!======================================================================
      Subroutine Read_astr(string,name,avalue)
!======================================================================
!     read characer value from string
!----------------------------------------------------------------------

      Implicit None

      Character(*) :: string,name, avalue
      Integer :: i,j

      i = INDEX(string,name); if(i.eq.0) Return 
      j = INDEX(string(i+1:),'=')
      if(j.gt.0) read(string(j+1:),*) avalue

      End Subroutine Read_astr

!======================================================================
      Subroutine Read_istr(string,name,ivalue)
!======================================================================
!     read integer value from string
!----------------------------------------------------------------------

      Implicit None

      Character(*) :: string,name
      Integer :: i,j, ivalue

      i = INDEX(string,name); if(i.eq.0) Return 
      j = INDEX(string(i+1:),'=')
      if(j.gt.0) read(string(j+1:),*) ivalue

      End Subroutine Read_istr

!======================================================================
      Subroutine Read_rstr(string,name,rvalue)
!======================================================================
!     read real value from string
!----------------------------------------------------------------------

      Implicit None

      Character(*) :: string,name
      Integer :: i,j
	  Real(8) :: rvalue

      i = INDEX(string,name); if(i.eq.0) Return 
      j = INDEX(string(i+1:),'=')
      if(j.gt.0) read(string(j+1:),*) rvalue

      End Subroutine Read_rstr
