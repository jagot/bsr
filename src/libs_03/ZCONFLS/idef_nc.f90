!======================================================================
      Integer Function Idef_ncfg(nu)
!======================================================================
!     gives the number of configuration in c-file (unit nu)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Integer :: ncfg
      Character(5) :: AS

      rewind(nu)
      ncfg=0
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      read(nu,'(a)') AS
      ncfg=ncfg+1
      go to 1
    2 rewind(nu)
      Idef_ncfg=ncfg

      End Function Idef_ncfg


!======================================================================
      Subroutine Jdef_ncfg(nu,ncfg,kcfg)
!======================================================================
!     gives the number of configuration in c-file (unit nu)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Integer, intent(out) :: ncfg,kcfg
      Character(5) :: AS

      rewind(nu)
      ncfg=0; kcfg = 0
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      ncfg=ncfg+1
      kcfg=kcfg+INDEX(AS,')',BACK=.TRUE.)/8
      read(nu,'(a)') AS
      go to 1
    2 rewind(nu)

      End Subroutine Jdef_ncfg
