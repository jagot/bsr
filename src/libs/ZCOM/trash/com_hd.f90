!======================================================================
      Subroutine Read_rarg(name,rvalue)
!======================================================================
!
!     read real argument
!
!----------------------------------------------------------------------

      Implicit None

      Character(*) :: name
      Real(8) :: rvalue

      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS
      Integer, External :: IARGC

      iarg = IARGC(); if(iarg.eq.0) Return
      iname=LEN_TRIM(name)
      Do i=1,iarg
       Call GETARG(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       read(AS(i1:i2),*) rvalue; Exit
      End do

      End Subroutine Read_rarg


!======================================================================
      Subroutine Read_iarg(name,ivalue)
!======================================================================
!
!     read integer argument
!
!----------------------------------------------------------------------

      Implicit None

      Character(*) :: name
      Integer :: ivalue

      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS
      Integer, External :: IARGC

      iarg = IARGC(); if(iarg.eq.0) Return
      iname=LEN_TRIM(name)
      Do i=1,iarg
       Call GETARG(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       read(AS(i1:i2),*) ivalue; Exit
      End do

      End Subroutine Read_iarg

!======================================================================
      Subroutine Read_aarg(name,avalue)
!======================================================================
!
!     read characer argument
!
!----------------------------------------------------------------------

      Implicit None

      Character(*) :: name, avalue

      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS
      Integer, External :: IARGC

      iarg = IARGC(); if(iarg.eq.0) Return
      iname=LEN_TRIM(name)
      Do i=1,iarg
       Call GETARG(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       read(AS(i1:i2),'(a)') avalue; Exit
      End do


      End Subroutine Read_aarg


!======================================================================
      Subroutine Read_iarr(name,na,iarr)
!======================================================================
!     read integer arrray from string name=a,b,c-d,f...
!----------------------------------------------------------------------

      Implicit None

      Character(*) :: name
      Integer :: na,iarr(na)

      Integer :: iarg,ia,iname,i,i1,i2,j,j1,j2,k,k1,k2
      Character(80) :: AS
      Integer, External :: IARGC

      iarg = IARGC(); if(iarg.eq.0) Return
      iname=LEN_TRIM(name)
      k=0
      Do i=1,iarg
       Call GETARG(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       k=1; Exit
      End do
      if(k.eq.0) Return

      ia=0; j1=i1; ! iarr=0
      Do
       j2=INDEX(AS(j1:i2),',')
       if(j2.eq.0) then; j2=i2; else; j2=j2+j1-1; end if
       j=0 ! INDEX(AS(j1:j2),'-'); k = j+j1-1
       if(j.eq.0) then
        ia=ia+1; if(ia.gt.na) Stop 'Read_iarr: ia > na'
        read(AS(j1:j2),*) iarr(ia)
       else
        read(AS(j1:k-1),*) k1
        read(AS(k+1:j2-1),*) k2
        Do k=k1,k2
         ia=ia+1; if(ia.gt.na) Stop 'Read_iarr: ia > na'
         iarr(ia)=k
        End do
       end if
       j1=j2+1;
       if(j1.gt.i2) Exit
      End do

      End Subroutine Read_iarr


!======================================================================
      Subroutine Read_name(name)
!======================================================================
!     read characer argument
!----------------------------------------------------------------------

      Implicit none

      Character(*) :: name

      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS
      Integer, External :: IARGC

      iarg = IARGC(); if(iarg.eq.0) Return
      Do i=1,iarg
       Call GETARG(i,AS)
       if(INDEX(AS,'=').ne.0) Cycle
       name=AS
       Exit
      End do

      End Subroutine Read_name
!======================================================================
      Subroutine Read_ipar(nu,name,ivalue)
!======================================================================
!
!     read the integer variable 'ivalue' with identifier 'name'
!     from unit 'nu', where the record like
!
!     name =  #####      is supposed to exist
!
!----------------------------------------------------------------------

      Implicit none

      Integer, Intent(in) :: nu
      Character(*), Intent(in) :: name
      Integer :: ivalue

      Character(80) :: AS
      Integer :: i

      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:i).ne.name) go to 1
      if(AS(i+1:i+1).ne.'='.and.AS(i+1:i+1).ne.' ') go to 1
      i=INDEX(AS,'=')+1
      read(AS(i:),*) ivalue
    2 Continue

      End Subroutine Read_ipar

!======================================================================
      Subroutine Read_rpar(nu,name,rvalue)
!======================================================================
!
!     read the real variable 'rvalue' with identifier 'name'
!     from unit 'nu', where the record like
!
!     name =  #####    ! coments
!
!     is supposed to exist
!
!----------------------------------------------------------------------

      Implicit none

      Integer, Intent(in) :: nu
      Character(*), Intent(in) :: name
      Real(8) :: rvalue

      Character(80) :: AS
      Integer :: i

      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:i).ne.name) go to 1
      if(AS(i+1:i+1).ne.'='.and.AS(i+1:i+1).ne.' ') go to 1
      i=INDEX(AS,'=')+1
      read(AS(i:),*) rvalue
    2 Continue

      End Subroutine Read_rpar


!======================================================================
      Subroutine Read_rval(nu,name,rvalue)
!======================================================================
!
!     read the real variable 'rvalue' with identifier 'name'
!     from unit 'nu', where the record like
!
!     name =  #####    ! coments
!
!     is supposed to exist
!
!----------------------------------------------------------------------

      Implicit none

      Integer, Intent(in) :: nu
      Character(*), Intent(in) :: name
      Real(8) :: rvalue

      Character(80) :: AS
      Integer :: i

      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      i=INDEX(AS,name)
      if(i.eq.0) go to 1
      i=INDEX(AS,'=')+1
      read(AS(i:),*) rvalue
    2 Continue

      End Subroutine Read_rval

!======================================================================
      Subroutine Read_apar(nu,name,avalue)
!======================================================================
!
!     read the character variable 'avalue' with identifier 'name'
!     from unit 'nu', where the record like
!
!     name =  #####    ! coments
!
!     is supposed to exist
!
!----------------------------------------------------------------------

      Implicit none

      Integer, Intent(in) :: nu
      Character(*), Intent(in) :: name
      Character(*) :: avalue

      Character(180) :: AS
      Integer :: i

      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:i).ne.name(1:i)) go to 1
!      if(AS(i+1:i+1).ne.'='.and.AS(i+1:i+1).ne.' ') go to 1
      i=INDEX(AS,'=')+1
      read(AS(i:180),*) avalue
    2 Continue

      End Subroutine Read_apar

!======================================================================
      Subroutine Read_iarray(nu,name,nv,ivalue)
!======================================================================
!
!     read the integer variable 'ivalue' with identifier 'name'
!     from unit 'nu', where the record like
!
!     name =  # # # # #      is supposed to exist
!
!----------------------------------------------------------------------

      Implicit none

      Integer,Intent(in) :: nu, nv
      Character(*),Intent(in) :: name
      Integer,Intent(out) :: ivalue(nv)

      Character(80) :: AS
      Integer :: i,j

      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      i = INDEX(AS,name)
      if(i.eq.0) go to 1
      i=INDEX(AS,'=')+1
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
!     find fisrt record in file containing "name"
!----------------------------------------------------------------------

      Implicit none

      Integer, Intent(in) :: nu
      Character(*), Intent(in) :: name

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
!     read the character variable 'avalue' with identifier 'name'
!     from unit 'nu', where the record like
!
!     name =  #####
!
!     is supposed to exist
!----------------------------------------------------------------------

      Implicit none

      Integer, Intent(in) :: nu
      Character(*), Intent(in) :: name
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
        !======================================================================
      Subroutine Check_file(AF)
!======================================================================

      Character(*), Intent(in) :: AF
      Logical :: EX

      Inquire(file=AF,exist=EX)
      if(.not.EX) then
       write(*,*) ' can not find file  ',AF;  Stop
      end if

      End Subroutine Check_file


!======================================================================
      Integer Function Icheck_file(AF)
!======================================================================

      Character(*), Intent(in) :: AF
      Logical :: EX

      Inquire(file=AF,exist=EX)
      Icheck_file = 1
      if(.not.EX) Icheck_file = 0

      End Function Icheck_file

!======================================================================
      Subroutine Find_free_unit(nu)
!======================================================================

      Implicit none
      Integer :: nu,i
      Logical :: connected

      nu = 0
      Do i=21,99
       INQUIRE(UNIT=i,OPENED=connected)
       if(connected) Cycle
       nu = i
       Exit
      End do
      if(nu.eq.0) Stop 'Find_free_unit: nu = 0'

      End Subroutine Find_free_unit
!--------------------------------------------------------------------
!     Electron-label block:
!
!     AL  - spectroscopic symbol for given L
!     LA  - value L from spectroscopic symbol
!     ELF4 - specroscopic notation for electron orbital (n,l,k) (a4)
!     ELF3 - specroscopic notation for electron orbital (n,l,k) (a3)
!     EL4_NLK(EL,n,l,k) -  decodes the specroscopic notation for a4
!     EL4_NLK(EL,n,l,k) -  decodes the specroscopic notation for a3
!--------------------------------------------------------------------


!====================================================================
      CHARACTER FUNCTION AL(L,k)
!====================================================================
!
!     provides spectroscopic symbols for L values
!
!     Limits: L <= 153
!
!--------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER(4), intent(in) :: L,K
      INTEGER(4) :: I
      CHARACTER(21), SAVE :: AS, AB

      DATA AS/'spdfghiklmnoqrtuvwxyz'/
      DATA AB/'SPDFGHIKLMNOQRTUVWXYZ'/

      i = L+1; IF(k.eq.5.or.k.eq.6) i=(L-1)/2+1
      if(i.ge.1.and.i.le.21) then
       if(k.eq.1.or.k.eq.5) AL=AS(I:I)
       if(k.eq.2.or.k.eq.6) AL=AB(I:I)
      elseif(i.ge.22.and.i.le.153) then
       AL=CHAR(i+101)  ! from {[} ...
      else
       write(*,*) 'L,k=',L,k
       i = i/0
       Stop ' AL: L out of range'
      end if

      END FUNCTION AL


!====================================================================
      Integer(4) Function LA(a)
!====================================================================
!     gives the value of L from spetroscopic symbol "a"
!--------------------------------------------------------------------

      IMPLICIT NONE
      Character, Intent(in) :: a
      CHARACTER(21), SAVE :: AS, AB
      Integer(4) :: i

      DATA AS/'spdfghiklmnoqrtuvwxyz'/
      DATA AB/'SPDFGHIKLMNOQRTUVWXYZ'/

      i = INDEX (AS,a)
      if(i.eq.0) i = INDEX (AB,a)
      if(i.eq.0) i = ICHAR(a)-101
      LA = i-1

      End Function LA


!====================================================================
      Character(4) function ELF4(n,l,k)
!====================================================================
!
!     gives the specroscopic notation for electron orbital (n,l,k)
!
!     n must be < 100; if they > 50, they are considered
!     as characters with code i+ICHAR(1)-1.
!
!     k must be < 61*61; if k<=61 - one character from ASET
!--------------------------------------------------------------------

      IMPLICIT NONE

      Integer(4), intent(in) :: n,l,k

      Integer(4) :: i,i1,k1,k2,kk

      Character(4) :: EL
      Character(1), EXTERNAL :: AL

      Character(61) :: ASET = &
	   '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

      EL='    '; i=4; i1=ICHAR('1')

      if(k.gt.0) then
       if(k.le.61) then
        EL(i:i)=ASET(k:k); i=i-1
       elseif(k.le.61*61) then
        k1=k/61; k2=mod(k,61);
        if(k2.eq.0) then; k1=k1-1; k2=61; end if
	       EL(i:i)=ASET(k2:k2); i=i-1
	       EL(i:i)=ASET(k1:k1); i=i-1
       else
        write(*,*) ' ELF4: k is out of limits:',k
        Stop
       end if
      end if

      EL(i:i)=AL(l,1);  i=i-1

      if(n.gt.0) then
       kk=n+i1-1
       if(kk.ge.65.and.kk.le.122) then           ! from A
        EL(i:i)=CHAR(kk)
       elseif(n.lt.10) then
        write(EL(i:i),'(i1)') n
       elseif(n.ge.10.and.n.lt.100.and.i.ge.2) then
        write(EL(i-1:i),'(i2)') n
       else
        EL(i:i) = 'n'
       end if
      end if

      ELF4=EL

      End Function ELF4


!====================================================================
      Character(3) function ELF3(n,l,k)
!====================================================================
!
!     gives the A3 specroscopic notation for electron orbital (n,l,k)
!
!     n and k must be < 10; if they > 9, they are replaced
!     by characters with code i+ICHAR(1)-1.
!
!--------------------------------------------------------------------

      IMPLICIT NONE

      Integer(4), intent(in) :: n,l,k

      Character(1), EXTERNAL :: AL

      Character(61) :: ASET = &
	   '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

      if(k.eq.0) then
       ELF3(1:1) = ' '
       ELF3(2:2)=CHAR(n + ICHAR('1') - 1)
       ELF3(3:3)=AL(l,1)
      else
       ELF3(1:1)=CHAR(n + ICHAR('1') - 1)
       ELF3(2:2)=AL(l,1)
       ELF3(3:3)=ASET(k:k)
      end if

      End Function ELF3


!====================================================================
      Subroutine EL4_NLK(EL,n,l,k)
!====================================================================
!
!     decodes the specroscopic notation for electron orbital (n,l,k)
!
!     It is allowed following notations: 1s, 2s3, 2p30, 20p3,
!     1sh, 20sh, kp, kp1, kp11, ns, ns3, ns33, ... .
!
!     Call:  LA, ICHAR
!--------------------------------------------------------------------
      IMPLICIT NONE

      Character(61) :: ASET = &
       '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

      Character(4), intent(in) :: EL
      Integer(4), Intent(out) :: n,l,k

      Integer(4) :: i,j,ic, k1,k2
      Integer(4), EXTERNAL :: LA

      n=0; l=-1; k=0; i=1; j=1
    1 if(EL(i:i).eq.' ') then
       i=i+1
       if(i.le.4) go to 1
      end if
      if(i.eq.4.and.j.eq.1) j=2
      if(i.gt.4) Return

      ic=ICHAR(EL(i:i))-ICHAR('1')+1

      if(j.eq.1) then                      !  n -> ?
       if(n.eq.0.and.ic.le.9) then
        n=ic; j=1
       elseif(n.eq.0.and.ic.gt.9) then
        n=ic; j=2
       elseif(ic.gt.9) then
        l=LA(EL(i:i)); j=3
       else
        n=n*10+ic
       end if

      elseif(j.eq.2) then                   !  l -> ?
       l=LA(EL(i:i)); j=3

      elseif(j.eq.3) then                   !  k -> ?
       if(i.eq.3) then
        k1 = INDEX(ASET,EL(i:i)); k2 = INDEX(ASET,EL(i+1:i+1))
        if(k2.gt.0) k = k1*61+k2
        if(k2.eq.0) k = k1
       else
        k = INDEX(ASET,EL(i:i))
       end if
       j=4
      end if

      i=i+1
      if(i.le.4.and.j.le.3) go to 1
      if(n.ge.100.or.l.lt.0) then
       write(*,*) 'EL4_nlk is fail to decode: ',EL
       Stop ' '
      end if

      End Subroutine EL4_NLK



!====================================================================
      Subroutine EL3_NLK(EL,n,l,k)
!====================================================================
!
!     decodes the specroscopic notation for electron orbital (n,l,k)
!
!     It is allowed following notations: 1s, 2s3, 20p3, 1sh, kp, kp1,
!                                        ns, ns3, nsa, ... .
!
!     Call:  LA, ICHAR
!--------------------------------------------------------------------

      IMPLICIT NONE

      Character(3), intent(in) :: EL
      Integer(4), intent(out) :: n,l,k

      Integer(4) :: i,j,ic
      Integer(4), EXTERNAL :: LA

      Character(61) :: ASET = &
	   '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

      n=0; l=-1; k=0; i=1
      j=1
    1 if(EL(i:i).eq.' ') then
       i=i+1
       if(i.le.3) go to 1
      end if
      if(i.gt.3) Return

      ic=ICHAR(EL(i:i))-ICHAR('1')+1

      if(j.eq.1) then            ! search n

       if(n.eq.0.and.ic.le.9) then
        n=ic
       elseif(n.eq.0.and.ic.gt.9) then
        n=ic
        j=2
       elseif(ic.gt.9) then
        l=LA(EL(i:i))
        j=3
       else
        n=n*10+ic
       end if

      elseif(j.eq.2) then        ! search l
       l=LA(EL(i:i));  j=3

      elseif(j.eq.3) then        ! search k

       k = INDEX(ASET,EL(i:i))

      end if

      i=i+1
      if(i.le.3.and.j.le.3) go to 1
      if(n.ge.100.or.l.lt.0.or.k.ge.62) &
       Stop 'EL4_NLK is fail to decode '

      End Subroutine EL3_NLK



!====================================================================
      Subroutine SORTA(n,S,IPT)
!--------------------------------------------------------------------
!     gives sorting pointer IPT for real(8) array S
!--------------------------------------------------------------------

      IMPLICIT none

      Integer(4), intent(in) :: n
      Real(8), intent(in), dimension(*) :: S
      Integer(4), intent(out), dimension(*) :: IPT

      Integer(4) :: i,i1,j1,i2,j2

      Do i=1,n; IPT(i)=i;  End do

      Do i1=1,n-1; j1=IPT(i1)
       Do i2=i1+1,n;  j2=IPT(i2)
        if(abs(S(j1)).lt.abs(S(j2))) then
         IPT(i2)=j1; j1=j2
        end if
       End do
       IPT(i1)=j1
      End do

      End Subroutine SORTA
