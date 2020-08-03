!======================================================================
      Subroutine Read_ipar(nu,name,ivalue)
!======================================================================
!     read integer variable 'ivalue' with identifier 'name'
!     from unit 'nu', where the record like    name = ...
!     is supposed to exist; if absent - ivalue is not changed 
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(*), intent(in) :: name
      Integer :: ivalue
      Character(80) :: AS
      Integer :: i

      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:i).ne.name) go to 1
      if(AS(i+1:i+1).ne.'='.and.AS(i+1:i+1).ne.' ') go to 1
      i=INDEX(AS,'=')+1
      if(i.eq.1) go to 1
      read(AS(i:),*) ivalue
    2 Continue

      End Subroutine Read_ipar 

    
!======================================================================
      Subroutine Read_rpar(nu,name,rvalue)
!======================================================================
!     read real variable 'rvalue' with identifier 'name'
!     from unit 'nu', where the record begining with  name = ...
!     is supposed to exist; if absent - rvalue is not changed 
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(*), intent(in) :: name
      Real(8) :: rvalue
      Character(80) :: AS
      Integer :: i

      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:i).ne.name) go to 1
      if(AS(i+1:i+1).ne.'='.and.AS(i+1:i+1).ne.' ') go to 1
      i=INDEX(AS,'=')+1
      if(i.eq.1) go to 1
      read(AS(i:),*,end=2,err=2) rvalue
    2 Continue

      End Subroutine Read_rpar 


!======================================================================
      Subroutine Read_rval(nu,name,rvalue)
!======================================================================
!     read real variable 'rvalue' with identifier 'name'
!     from unit 'nu', where the record containing  name = ...
!     is supposed to exist; if absent - rvalue is not changed 
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(*), intent(in) :: name
      Real(8) :: rvalue
      Character(80) :: AS
      Integer :: i

      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      i=INDEX(AS,name)  
      if(i.eq.0) go to 1
      i=INDEX(AS,'=')+1
      if(i.eq.1) go to 1
      read(AS(i:),*) rvalue
    2 Continue

      End Subroutine Read_rval 


!======================================================================
      Subroutine Read_apar(nu,name,avalue)
!======================================================================
!     read character variable 'avalue' with identifier 'name'
!     from unit 'nu', where the record begining with  name = ...
!     is supposed to exist; if absent - avalue is not changed 
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(*), intent(in) :: name
      Character(*) :: avalue
      Character(180) :: AS
      Integer :: i

      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:i).ne.name(1:i)) go to 1
      if(AS(i+1:i+1).ne.'='.and.AS(i+1:i+1).ne.' ') go to 1
      i=INDEX(AS,'=')+1
      if(i.eq.1) go to 1
      read(AS(i:180),*) avalue
    2 Continue

      End Subroutine Read_apar 
        

!======================================================================
      Subroutine Read_iarray(nu,name,nv,ivalue)
!======================================================================
!     read integer array 'ivalue(1:nv)' with identifier 'name'
!     from unit 'nu', where the record like    name = ...
!     is supposed to exist; if absent - ivalue is not changed 
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu, nv
      Character(*), intent(in) :: name
      Integer, intent(out) :: ivalue(nv)
      Character(80) :: AS
      Integer :: i,j
 
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      i = INDEX(AS,name)
      if(i.eq.0) go to 1
      i=INDEX(AS,'=')+1
      if(i.eq.1) go to 1
      read(AS(i:),*) (ivalue(j),j=1,nv)
    2 Continue

      End Subroutine Read_iarray 


!======================================================================
      Integer Function Ifind_position(nu,name)
!======================================================================
!     find position of line with "name" in the begining
!----------------------------------------------------------------------
      Implicit none
      Integer, Intent(in) :: nu
      Character(*), Intent(in) :: name
      Character(80) :: AS
      Integer :: i,j
 
      Ifind_position = 0
      i=LEN_TRIM(name); j=0
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      j=j+1
      if(AS(1:i).ne.name) go to 1
      Ifind_position = j
      Backspace(nu)
      Return
    2 rewind(nu)

      End Function Ifind_position


!======================================================================
      Integer Function Jfind_position(nu,name)
!======================================================================
!     find fisrt record in file containing "name"    not-finished yet?
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(*), intent(in) :: name
      Character(80) :: AS
      Integer :: i,j
 
      Jfind_position = 0
      i=LEN_TRIM(name); j=0
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      j=j+1
      if(AS(1:i).ne.name) go to 1
    2 Jfind_position = j
      Backspace(nu)

      End Function Jfind_position


!======================================================================
      Subroutine Read_string(nu,name,avalue)
!======================================================================
!     read character variable 'avalue' with identifier 'name'
!     from unit 'nu', where the record begining with  name = ...
!     is supposed to exist; if absent - avalue is not changed. 
!     Avalue may contains blanks.
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(*), intent(in) :: name
      Character(*) :: avalue
      Character(180) :: AS
      Integer :: i,j

      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:i).ne.name(1:i)) go to 1
      i=INDEX(AS,'=')+1
      j=LEN_TRIM(AS)
      avalue = AS(i:j)
    2 Continue

      End Subroutine Read_string
        