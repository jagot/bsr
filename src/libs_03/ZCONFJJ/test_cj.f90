!====================================================================
      Subroutine test_cj
!====================================================================
!     check the configuration in jj-coupling
!--------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer :: i,k, JW,JQ
      Integer, external :: Jterm, ITRA

! ... check the total number of electron:

      k = 0; Do i=1,no; k=k+iq(i); End do
      if(ne.eq.0) ne = k
      if(ne.ne.k) then
       write(*,*) 'Incorrect number of electrons, ne =',k,ne
       Call Incode_confj
       write(*,*) CONFIG(1:ia)
       Stop 'Stop in TEST_cj'
      end if

! ... check the total parity:

      k=0; Do i=1,no; k=k+iq(i)*ln(i);  End do
      k=(-1)**k
      if(parity.eq.0) parity = k
      if(k.ne.parity) then
       write(*,*) 'Incorrect total parity =',k,parity
       Call Incode_confj
       write(*,*) CONFIG(1:ia)
       Stop 'Stop in TEST_cj'
      end if

! ... check the number of electrons in the shell ...

      Do i=1,no
       if(iq(i).gt.jn(i)+1) then
        write(*,*) 'Shell ',i,' contains too many electrons'
        Call Incode_confj
        ia = LEN_TRIM(CONFIG)
        write(*,*) CONFIG(1:ia)
        Stop 'Stop in TEST_cj'
       end if
       if(nn(i).le.0) Stop ' TEST_cj: n-value < 0'
      End do

! ... check that subshell term is an allowed one ...

      Do i=1,no
       k=Jterm(jn(i),iq(i),0,Jshell(i),Vshell(i),JW,JQ)
       if(k.eq.0) then
        write(*,*) 'Term of shell ',i,' is not recognized'
        Call Incode_cj
        write(*,*) CONFIG(1:ia)
        write(*,*) SHELLJ(1:ia)
        Stop 'Stop in TEST_cj'
       end if
      End do

! ... check the coupling scheme:

       Do i=2,no
        k=ITRA(Jintra(i-1),Jshell(i),Jintra(i))
        if(k.eq.0) then
         write(*,*) 'The coupling of shell',i,' is incorrect'
         Call Incode_cj
         write(*,'(a)') CONFIG(1:ia)
         write(*,'(a)') SHELLJ(1:ia)
         write(*,'(9x,a)') INTRAJ(1:ia)
         Stop 'Stop in TEST_cj'
        end if
       End do

      End Subroutine test_cj


!======================================================================
      Subroutine test_ac
!======================================================================
!     check the "AFTER" relations between orbitals in configurations
!----------------------------------------------------------------------
      Use conf_jj; Use orb_jj

      Implicit none
      Integer :: ic, i,j, i1,i2
      Integer, external :: Ifind_jjorb, after

      Call Alloc_after_conditions(0)
      Do ic = 1,ncfg
       Call Get_cfg_jj(ic)

       Do i=1,no-1;   i1 = Ifind_jjorb(nn(i),kn(i),in(i),1)
       Do j=i+1,no;   i2 = Ifind_jjorb(nn(j),kn(j),in(j),1)

       if(i1.eq.i2) then
        write(*,'(/a)')    ' The same electron are on different shells:'
        write(*,'(/a,i6)') ' TEST_ac: configuration ', ic
        Call Print_conf_jj(0,0,0.d0)
        write(*,'(/a,a)')  ' Orbital: ', ELF(i1)
        Stop
       end if
       
       if(after(i1,i2).lt.0) then
        write(*,'(/a,i6)') ' TEST_ac: configuration ', ic
        Call Print_conf_jj(0,0,0.d0)
        write(*,'(/a,2a5)')  &
        ' The AFTER relations are not satisfied for orbitals: ', &
          ELF(i1), ELF(i2)
       end if

       End do
       End do

      End do   ! over ic

      End Subroutine test_ac



!======================================================================
      Subroutine test_aj
!======================================================================
!     check the "AFTER" relations between orbitals 
!----------------------------------------------------------------------
      Use conf_jj; Use orb_jj

      Implicit none
      Integer :: i,j, i1,i2
      Integer, external :: Ifind_jjorb, after

      Do i=1,no-1;   i1 = Ifind_jjorb(nn(i),kn(i),in(i),1)
      Do j=i+1,no;   i2 = Ifind_jjorb(nn(j),kn(j),in(j),1)

      if(i1.eq.i2) then
       write(*,'(/a)')    ' The same electron are on different shells:'
       Call Print_conf_jj(0,0,0.d0)
       write(*,'(/a,a)')  ' Orbital: ', ELF(i1)
       Stop
      end if
      
      if(after(i1,i2).lt.0) then
       Call Print_conf_jj(0,0,0.d0)
       write(*,'(/a,2a5)')  &
       ' The AFTER relations are not satisfied for orbitals: ', &
         ELF(i1), ELF(i2)
      end if

      End do
      End do

      End Subroutine test_aj




