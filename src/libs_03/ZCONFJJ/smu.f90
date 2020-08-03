!=======================================================================
      Real(8) Function smu (ka,kb,kc,kd,L,v,S)
!=======================================================================
!     Computes the coefficient for Breit operator, according to
!     Grant and Pyper, J.Phys.B 9, 761 (1976)
!     Grant book: p.345
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: ka,kb,kc,kd,L,v
      Real(8), intent(out) :: S(8)
      Integer :: i, la,lb,lc,ld, K,Kp, vv
      Integer, external :: ITRA, l_kappa
      Real(8) :: b,c,bp,cp

      SMU = 0.d0; S = 0.d0; if(L.lt.0) Return;  if(v.lt.0) Return 

      la=l_kappa(ka); lc=l_kappa(kc)
      if(mod(la+lc+v,2).eq.0) Return

      lb=l_kappa(kb); ld=l_kappa(kd)
      if(mod(lb+ld+v,2).eq.0) Return

      if(v.eq.L) then
       if(L.eq.0) Return
       SMU = -(ka+kc)*(kb+kd); SMU=SMU/L/(L+1); S=SMU
      else
       K=kc-ka; Kp=kd-kb
       if(v.eq.L-1) then
        bp=L+1; bp=bp/2/(L+L-1); cp=-(L-2); cp=cp/2/L/(L+L-1)
        S(1) = (L+K )*(bp+cp*KP)
        S(2) = (L+KP)*(bp+cp*K )
        S(3) = (L-K )*(bp-cp*KP)
        S(4) = (L-KP)*(bp-cp*K )
        S(5) =-(L+K )*(bp-cp*KP)
        S(6) =-(L-KP)*(bp+cp*K )
        S(7) =-(L-K )*(bp+cp*KP)
        S(8) =-(L+KP)*(bp-cp*K )
       elseif(v.eq.L+1) then
        i = -L-1
        b=L; b=b/2/(L+L+3); c=(L+3); c=c/(L+L+2)/(L+L+3)
        S(1) = (i+KP)*(b+c*K )
        S(2) = (i+K )*(b+c*KP)
        S(3) = (i-KP)*(b-c*K )
        S(4) = (i-K )*(b-c*KP)
        S(5) =-(i-KP)*(b+c*K )
        S(6) =-(i+K )*(b-c*KP)
        S(7) =-(i+KP)*(b-c*K )
        S(8) =-(i-K )*(b+c*KP)
       end if

       Do i=1,8; SMU = SMU + abs(S(i)); End do

      end if

      End Function smu 


