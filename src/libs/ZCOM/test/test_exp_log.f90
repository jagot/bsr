    Implicit real(8) (A-H,O-Z)
    t = 1.d0;  a = 1.d0
    write(*,'(a,f10.6,5x,a,2f10.6)') 't=',t,'a=',a,EXP_LOG(a,t)
    write(*,'(a,f10.6,5x,a,2f10.6)') 't=',t,'a=',a,EXP_LOG(a,t) * exp(1.d0/t)
    END
