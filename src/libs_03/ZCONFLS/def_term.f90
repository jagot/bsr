!====================================================================
      Subroutine Def_term(nu,ILT,IST,IPT)
!====================================================================
!     define the term from c-file (unit nu) in (2j+1) format
!--------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer, intent(in) :: nu
      Integer, intent(out) :: ILT,IST,IPT
      Integer :: i
      Character(80) :: AS

! ... LS case:

      rewind(nu)
    1 Read(nu,'(a)') CONFIG; if(CONFIG(5:5).ne.'(') go to 1
      Read(nu,'(a)') COUPLE
      Call DECODE_c
      IPT=SUM(ln(1:no)*iq(1:no)); IPT=(-1)**IPT
      ILT=LS(no,4); IST=LS(no,5) 

! ... check if it is a J-case:

      rewind(nu)
      read(nu,'(a60)') AS;  i=INDEX(AS,'=')
      if(i.gt.0) then
       read(AS(i+1:),*) ILT; IST=0
      end if

      End Subroutine DEF_term


!====================================================================
      Subroutine Def_term_BSR(nu,ILT,IST,IPT)
!====================================================================
!     define the term from c-file (unit nu) as L,(2J), 2S+1
!--------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer, intent(in) :: nu
      Integer, intent(out) :: ILT,IST,IPT
      Integer :: i
      Character(80) :: AS

! ... LS case:

      rewind(nu)
    1 Read(nu,'(a)') CONFIG; if(CONFIG(5:5).ne.'(') go to 1
      Read(nu,'(a)') COUPLE
      Call DECODE_c
      IPT=SUM(ln(1:no)*iq(1:no)); IPT=(-1)**IPT
      ILT=(LS(no,4)-1)/2; IST=LS(no,5) 

! ... check if it is a J-case:

      rewind(nu)
      read(nu,'(a60)') AS;  i=INDEX(AS,'=')
      if(i.gt.0) then
       read(AS(i+1:),*) ILT; IST=0
      end if

      End Subroutine DEF_term_BSR
