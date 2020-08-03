!====================================================================
      Subroutine read_peel_orb(nu)
!====================================================================
!     reads the common closed shells list from GRASP c-file (unit nu)
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(280) :: AS
      Character(5) :: EL
      Integer :: k, n,l,j,iset,ka
      Integer, external :: Ifind_jjorb

      rewind(nu)
      Do 
       read(nu,'(a)',end=1) AS
       if(AS(1:1).eq.'*') Exit
       if(AS.ne.'Peel subshells:') Cycle
       Do 
        read(nu,'(a)',end=1) AS
        if(AS(1:1).eq.'*') Exit
        if(AS.eq.'CSF(s):') go to 1
        k=1
        Do 
         EL = AS(k:k+4); if(len_trim(EL).eq.0) Exit
         Call EL_NLJK(EL,n,ka,l,j,iset)
         j=Ifind_jjorb(n,ka,iset,2)
         k = k + 5
        End do
       End do
      End do

    1 Continue

      End Subroutine read_peel_orb

