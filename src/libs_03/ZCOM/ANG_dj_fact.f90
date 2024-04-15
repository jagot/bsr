!======================================================================
      Real(8) Function DJ_fact(IL1,IL2,IS1,IS2,JOT1,JOT2,k)
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
      Integer, intent(in) :: IL1,IL2,IS1,IS2,JOT1,JOT2,k
      Real(8), external :: Z_6j

      DJ_fact = 0.d0
      if(IS1.eq.IS2) DJ_fact = (-1)**((IS1+IL1+JOT2-3)/2+k)*  &
                               Dsqrt(1.d0*JOT1*JOT2)*         &
                               Z_6j(IL1,JOT1,IS1,JOT2,IL2,k+k+1)

      End Function DJ_fact


!====================================================================
      Real(8) Function DJM_fact(IL1,IL2,IS1,IS2,JOT1,JOT2,k)
!====================================================================
!     defines the J-dependence for reduced matrix elements of
!     k-pole magnetic operator MB:
!
!     <LSJ||MB[k]||L'S'J'> = DJM_fact *  <LS||MB[k]||L'S'>
!
!     DJM_fact =  [J,J',k]^1/2  { L  S  J }
!                               { L' S' J'}
!                               {k-1 1  k }
!
!     momenta are given in (2J+1)-representation
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: IL1,IL2,IS1,IS2,JOT1,JOT2,k
      Real(8), external :: Z_9j

      DJM_fact = Dsqrt(1.d0*JOT1*JOT2*(k+k+1))*         &
                 Z_9j(IL1,IS1,JOT1,IL2,IS2,JOT2,k+k-1,3,k+k+1)

      End Function DJM_fact
