!======================================================================
!     defines the J-dependence for reduced matrix elements of
!     k-pole electric and magnetic transition operator (type MA):
!
!     <LSJ||O[k]||L'S'J'> = DJ_fact *  <LS||O[k]||L'S'>  
!
!     DJ_fact = (-1)^(L+S+J'+k) [J,J']^1/2  { L  S J },  S = S' 
!                                           { J' k L'}
!
!     momenta are given in (2J+1)-representation
!======================================================================


      Implicit none
      Integer :: IL1,IL2,IS1,IS2,JOT1,JOT2,k, inp,out
      Real(8) :: S
      Real(8), external :: Z_6j

      inp=1;  open(inp, file = 'test_LSJ.inp')
      out=2;  open(out, file = 'test_LSJ.out')


    1 Continue
      read(inp,*,end=2) k,IS1,IL1,JOT1,IS2,IL2,JOT2 
      if(k.lt.0) Stop
      S = Z_6j(2*IL1+1,JOT1+1,IS1,JOT2+1,2*IL2+1,k+k+1) 
      S = S*S * (2*IL1+1)*(JOT2+1)
      if(IS1.ne.IS2) S = 0.d0      
      write(out,'(7i5,f10.5)') k,IS1,IL1,JOT1,IS2,IL2,JOT2, S
      go to 1
    2 Continue

      End
