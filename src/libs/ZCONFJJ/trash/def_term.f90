!====================================================================
      Subroutine Def_JP(nu,J,P)
!====================================================================
!
!     define the term from GRASP c-file (unit nu)
!
!--------------------------------------------------------------------

      USE conf_jj

      IMPLICIT NONE

      Integer(4), INTENT(in) :: nu
      Integer(4), INTENT(out) :: J,P

      J=0; P=0

! ... GRASP case:

      rewind(nu)
      read(nu,'(a)') AS
      if(AS(1:4).eq.'Core') then
    1  Read(nu,'(a)') CONFIG; if(CONFIG(6:6).ne.'(') go to 1
       Read(nu,'(a)') SHELLJ
       Read(nu,'(3x,a)') INTRAJ
       Call Decode_cj
       P=SUM(ln(1:no)*iq(1:no)); P=(-1)**P
       J=Jintra(no)
      end if

      End Subroutine DEF_JP

