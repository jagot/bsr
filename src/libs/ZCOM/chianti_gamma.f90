!------------------------------------------------------------
      Real(8) Function Chianti_gamma(T,n,x,y,met,C,de)
!------------------------------------------------------------ 
!     Gamma-value from Chianti format at temperature T
!     T - temperature in K
!     de - transitions energy (in Ry)
!------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,met
      Real(8), intent(in) :: T,C,de,x(n),y(n)
      Real(8), parameter :: Ry = 13.6057, K_eV = 11604.5 
      Real(8) :: e,xe,xx,yy
      Real(8) , external :: XLAGR

      xe = T/K_eV/Ry/de
      e = exp(1.d0)

      Select Case(met)
       Case(1,4);     xx = 1.d0 - log(C)/log(xe+C)
       Case(2,3);     xx = xe/(xe+C)
      End Select

      yy = XLAGR (n,n,x,y,xx)

      Select Case(met)
       Case(1)
        Chianti_gamma = yy * log(xe+e)
       Case(2) 
        Chianti_gamma = yy 
       Case(3)
        Chianti_gamma = yy / (xe + 1.d0)
       Case(4)
        Chianti_gamma = yy * log(xe+C)
      End Select
        
      End Function Chianti_gamma 

      
                                      