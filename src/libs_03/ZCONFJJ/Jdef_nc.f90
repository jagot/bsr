!======================================================================
      Integer Function Jdef_ncfg(nu)
!======================================================================
!     defines the number of configuration in c-file (unit nu)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Integer :: ncfg
      Character(6) :: AS

      rewind(nu)
      ncfg=0
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(6:6).ne.'(') go to 1
      ncfg=ncfg+1
      go to 1
    2 rewind(nu)
      Jdef_ncfg=ncfg

      End Function Jdef_ncfg


!======================================================================
      Subroutine Jdef_kcfg(nu,ncfg,kcfg)
!======================================================================
!     defines the number of configuration in c-file (unit nu),
!     plus the total "length" of these configurations, kcfg
!----------------------------------------------------------------------
      Implicit none
      Integer :: i,nu,ncfg,kcfg
      Character(300) :: AS

      rewind(nu)
      ncfg=0; kcfg=0
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(6:6).ne.'(') go to 1
      ncfg=ncfg+1
      i = INDEX(AS,')',BACK=.TRUE.)
      if((i/9)*9.ne.i) Stop 'Jdef_ncfg: problems with kcfg'
      kcfg = kcfg + i/9
      go to 1
    2 rewind(nu)

      End Subroutine Jdef_kcfg


