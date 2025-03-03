!======================================================================
      Subroutine Read_conf_jj
!======================================================================
!     Read the configuration list from c-file (unit 'nuc'),
!     define input angular symmetries and compare with existing ones
!     (in file with unit 'nub')
!----------------------------------------------------------------------

      USE param_jj

      USE conf_jj
      USE symc_list
      USE symt_list

      Implicit none

      Integer, External :: Jdef_ncfg, Iadd_symc, Iadd_symt, &
                           Iort_conf_jj, DEF_ij
      Integer :: i,j,k,ic,jc,ik,jk,it,jt,ii,ij,ijc, iort_c

! ... define ncfg:

      Call Alloc_cfg (0)
      Call Alloc_symc(0)
      Call Alloc_symt(0)

      ncfg=Jdef_ncfg(nuc);
      if(ncfg.eq.0) Stop ' ncfg = 0, nothing to do '

! ... read old information, if any:

      if(new) then
       Call Alloc_symc(0)
       Call Alloc_symt(0)
      else
       Call Read_symc(nub)
       Call Read_symt(nub)
!       IT_conf = - IT_conf                ???
      end if

!---------------------------------------------------------------------
! ... define new symmetries from c-file:

      rewind(nuc); ne=0; parity=0
    1 read(nuc,'(a)',end=2) CONFIG
      if(CONFIG(1:3).eq.'*') go to 2
      if(CONFIG(6:6).ne.'(') go to 1
      read(nuc,'(a)') SHELLJ
      read(nuc,'(5x,a)') INTRAJ

      Call Decode_cj; Call Test_cj

! ... check angular symmetry:

      iconf = Iadd_symc(Jtotal,no,iq,kn)
      iterm = Iadd_symt(iconf,no,Jshell,Vshell,Jintra)

      go to 1
    2 Continue     

      write(pri,'(a,i5)') ' number of atomic states   = ',NCFG
      write(pri,'(a,i5)') ' number of configurations  = ',NSYMC
      write(pri,'(a,i5)') ' number of ang. symmetries = ',NSYMT
      write(pri,*)

!----------------------------------------------------------------------
! ... define IP_term and IC_term pointers:

      if(Allocated(IT_stat)) Deallocate(IT_stat)
                             Allocate(IT_stat(nsymt))
      if(Allocated(JP_term)) Deallocate(JP_term)
                             Allocate(JP_term(nsymt))
      if(Allocated(IC_term1)) Deallocate(IC_term1)
                              Allocate(IC_term1(nsymc))
      if(Allocated(IC_term2)) Deallocate(IC_term2)
                              Allocate(IC_term2(nsymc))

      IT_stat = 1
      Do i=1,nsymt
       if(IT_conf(i).gt.0) Cycle
       IT_conf(i) = iabs(IT_conf(i))
       IT_stat(i) = 0
      End do

      Call SortI(nsymt,IT_conf,JP_term)

      IC_term1=0
      Do i=1,nsymt
       ic=IT_conf(JP_term(i))
       if(IC_term1(ic).eq.0) IC_term1(ic)=i
                             IC_term2(ic)=i
      End do

!----------------------------------------------------------------------
! ... define if we need additional calculations

      if(allocated(IC_need)) Deallocate(IC_need)
                             Allocate(IC_need(nsymc))
      if(allocated(JC_need)) Deallocate(JC_need)
                             Allocate(JC_need(nsymc*(nsymc+1)/2))
      if(allocated(IT_need)) Deallocate(IT_need)
                             Allocate(IT_need(nsymt))
      if(allocated(IT_done)) Deallocate(IT_done)
                             Allocate(IT_done(nsymt*(nsymt+1)/2))

      IT_done=0
      if(new) then
        icalc=.TRUE.; IC_need=1; JC_need=1; IT_need=1
      else
       Call Read_done(nub)

       icalc=.FALSE.; IC_need = 0; JC_need = 0; IT_need = 0
       Do ic = 1,nsymc
        Call Get_symc(ic,Jtotal1,no1,nn1,kn1,ln1,jn1,iq1,in1) 
       Do jc = 1,ic
        Call Get_symc(jc,Jtotal2,no2,nn2,kn2,ln2,jn2,iq2,in2)      
         iort_c = Iort_conf_jj(2) 
         if(Jtotal1.ne.Jtotal2) iort_c = 1

         k = 0 
         Do ik=IC_term1(ic),IC_term2(ic);  it=JP_term(ik)
         Do jk=IC_term1(jc),IC_term2(jc);  jt=JP_term(jk)
           ij = DEF_ij(it,jt)
           if(IT_done(ij).eq.0) IT_done(ij) = iort_c
           if(IT_stat(it)*IT_stat(jt).eq.0.and.IT_done(ij).eq.0) &
              IT_done(ij)=-1
           if(IT_done(ij).eq.0) k = 1
         End do
         End do
         if(k.eq.0) Cycle
         ijc=DEF_ij(ic,jc); JC_need(ijc)=1; IC_need(ic)=1; IC_need(jc)=1
         icalc=.TRUE.
       End do
       End do

       ! ... define IT_need:

       Do it = 1,nsymt; k = 0
        Do jt = 1,nsymt
         ij=DEF_ij(it,jt); if(IT_done(ij).ne.0) Cycle; k=1; Exit
        End do         
        IT_need(it) = k
       End do

      end if

      if(icalc) then
       write(*,'(/a,a4/)')   'Need of additional calculations --> yes '
      else
       write(*,'(/a,a4/)')   'Need of additional calculations --> no '
      end if
      
      Deallocate(IT_stat)

      End Subroutine read_conf_jj


