!--------------------------------------------------------------------
      Integer(4) FUNCTION IWF(NN,LL,KK)
!--------------------------------------------------------------------
!
!     find the index for orbital (n,l,k) in RADIAL list,
!     otherwise stop
!
!--------------------------------------------------------------------

      USE RADIAL

      Integer(4), Intent(in) :: NN,LL,KK

      Do i=1,NRF
       if( NN.eq.nro(i).and.LL.eq.lro(i).and.KK.eq.kro(i) ) then
        IWF=i;  Return
       end if
      End do

      write(*,'(/'' IWF --> cannot find the orbital:'', &    
         5x,''N='',i3,3x,''L='',i3,3x,''K='',i3)')  NN,LL,KK  
      Stop

      END FUNCTION IWF



!--------------------------------------------------------------------
      Integer(4) FUNCTION JWF(NN,LL,KK)
!--------------------------------------------------------------------
!
!     find the index for orbital (n,l,k) in RADIAL list,
!     otherwise returns 0
!
!--------------------------------------------------------------------

      USE RADIAL

      Integer(4), Intent(in) :: NN,LL,KK

      JWF=0
      Do i=1,NRF
       if( NN.eq.nro(i).and.LL.eq.lro(i).and.KK.eq.kro(i) ) then
        JWF=i; Return
       end if
      End do

      END FUNCTION JWF
