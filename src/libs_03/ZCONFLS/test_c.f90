!====================================================================
      Subroutine TEST_c
!====================================================================
!     check the configuration
!--------------------------------------------------------------------
      Use conf_LS, iiterm => iterm

! ... check the total number of electron ...

      k = 0
      Do i=1,no; k=k+iq(i);  End do
      if(ne.eq.0) ne = k
      if(ne.ne.k) then
       write(*,*) ' Incorrect number of electrons !',k,ne
       Call Pri_conf(0,0,0.d0)
       Stop
      end if

! ... check the total parity ...

      k=0
      Do i=1,no
       k=k+iq(i)*ln(i)
      End do
       k=(-1)**k
       if(parity.eq.0) parity = k
       if(k.ne.parity) then
        write(*,'(/'' Incorrect total parity !''/)')
        Call Pri_conf(0,0,0.d0)
        Stop
       end if

! ... check the number of electrons in the shell ...

      Do i=1,no
       if(iq(i).gt.4*ln(i)+2) then
        write(*,'(/a,i2,a)') ' Shell ',i,' contains too many electrons'
        Call Pri_conf(0,0,0.d0)
        Stop
       end if
       if(nn(i).le.0) Stop ' TEST_C: n-value < 0'
      End do

! ... check that subshell term is an allowed one ...

      Do i=1,no
       k = Iterm_LS(ln(i),iq(i),0,LS(i,1),LS(i,2),LS(i,3))
       if(k.eq.0) then
        write(*,'(/a,i2,a)') ' Term of shell ',i,' is not recognized'
        Call Pri_conf(0,0,0.d0)
        Stop
       end if
      End do

! ... check the coupling scheme ...

       Do i=2,no
        kl = ITRI(LS(i-1,4),LS(i,2),LS(i,4))
        ks=ITRI(LS(i-1,5),LS(i,3),LS(i,5))
        if(kl*ks.eq.0) then
         write(*,'(/a,i2,a)') ' The coupling of shell',i,' is incorrect'
         Call Pri_conf(0,0,0.d0)
         Stop
        end if
       End do

      End Subroutine TEST_c
