!======================================================================
      Subroutine TEST_a
!======================================================================
!     check the AFTER relations between orbitals for all configurations
!----------------------------------------------------------------------
      Use conf_LS; Use orb_LS  

      Do ic = 1,ncfg
      
      Call Get_cfg_LS(ic)
      ip = IP_state(ic)

       Do i=1,no-1; i1 = IP_orb(ip+i)
       Do j=i+1,no; i2 = IP_orb(ip+j)

       if(i1.eq.i2) then

        IORT(i1,i2) = 1

       else

       j1=min(i1,i2)
       j2=max(i1,i2)
       jp=(i2-i1)/IABS(i2-i1)
       if(IORT(j1,j2).eq.0) then
        IORT(j1,j2)=jp
       elseif(IORT(j1,j2).ne.jp) then
        write(*,'(/a,i6)') ' TEST_a: configuration ', ic
        Call Pri_conf(0,0,0.d0)
        write(*,'(/a,2a5,i5)')  &
          ' The AFTER relations are not satisfied for orbitals: ', &
            ELF(i1), ELF(i2), IORT(j1,j2)
       end if

       end if

       End do
       End do

      End do   ! over ic

      End Subroutine TEST_a


