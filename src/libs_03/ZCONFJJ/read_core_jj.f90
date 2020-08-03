!====================================================================
      Subroutine read_core_jj(nu)
!====================================================================
!     reads the common closed shells list from GRASP c-file (unit nu)
!--------------------------------------------------------------------
      Use conf_jj
      
      Implicit none
      Integer, intent(in) :: nu
      Character(5) :: EL
      Integer, external :: Ifind_jjorb
      Integer :: i,ii, n,l,j,iset,ka

      rewind(nu)
      read(nu,'(/a)') core
      ncore=0
      ii=1
      Do i=1,mcore
       EL=core(ii:ii+4);  if(EL(3:5).eq.'  ') Exit
       ncore=ncore+1
       Call EL_NLJK(EL,n,ka,l,j,iset)
       nn_core(ncore)=n; l_core(ncore)=l; j_core(ncore)=j
       k_core(ncore) = ka
       if(iset.ne.0) Stop 'R_core: iset > 0'
       j=Ifind_jjorb(n,ka,iset,2)
       ii=ii+5
      End do

      End Subroutine read_core_jj

