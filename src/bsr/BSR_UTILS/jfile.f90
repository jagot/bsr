!======================================================================
!     UTILITY  J F I L E
!
!                   C O P Y R I G H T -- 2007
!
!     Written by:   Oleg Zatsarinny
!
!----------------------------------------------------------------------
!     find the maximum contributions for given configurations
!     (first c-file) from the list of solutions from j-file
!     Results in the name.cc file
!
!     Example:  jfile name=1 sol=1-5 eps_c=0.0001
!
!              - analysing the first 5 solutions in
!                1.j and output in 1.c all configurations
!                in order of importance
!                Now you may delete all configurations with
!                small coefficients and then should (!)
!                run again CI program because order in c-file
!                does not agree with old l-file
!----------------------------------------------------------------------
      Use conf_jj

      Implicit real(8) (A-H,O-Z)
      Character(80) :: name = '?', sol = '?' , AF

      Integer :: nuc = 1       !   name.c
      Integer :: nuj = 2       !   name.j
      Integer :: nur = 3       !   name.cc

      Real(8) :: eps_c = -1.d0
      Real(8) :: cc = 1.d0

      Integer, allocatable :: IP(:), IS(:)
      Real(8), allocatable :: C(:)
      Integer, external :: Ipointer, Ifind_position

!----------------------------------------------------------------------
! ... input data:

      Call Read_aarg('name',name)
      Call Read_aarg('sol',sol)
      if(name(1:1).eq.'?'.or.sol(1:1).eq.'?') &
       Stop 'Call as  jfile name=... sol=1,3,4-5... eps_c=...'
      Call Read_rarg('eps_c',eps_c)
      Call Read_rarg('cc',cc)

! ... input configurations:

      iname=LEN_TRIM(name)

      AF=name(1:iname)//'.c'; Call Check_file(AF);  Open(nuc,file=AF)

      Call Read_conf_jj(nuc,0,'add','nocheck')

      Allocate(C(ncfg),IP(ncfg),IS(ncfg));   C=0.d0;  WC=0.d0; IS=0

      Call Read_iarr('sol',ncfg,IS)

! ... proceed solutions:

      AF=name(1:iname)//'.j'; Call Check_file(AF);  Open(nuj,file=AF)

      Call Read_ipar(nuj,'ncfg',i)
      if(i.ne.ncfg) Stop 'jfile:  ncfg in j-file <> ncfg in c-file'
      Call Read_ipar(nuj,'nsol',nsol)

      i=Ifind_position(nuj,'Solutions');  read(nuj,*)

      Do i=1,nsol
       read(nuj,*) j
       read(nuj,*) E,jj,ic1,ic2
       read(nuj,*) C(ic1:ic2)
       if(Ipointer(ncfg,IS,j).eq.0) Cycle
       Do ic=ic1,ic2
        S = abs(C(ic)); if(S.gt.WC(ic)) WC(ic)=S
       End do
      End do

! ... write the results in ordered form:


      AF=name(1:iname)//'.cc';  Open(nur,file=AF)

      rewind(nuc)
      Do
       read(nuc,'(a)') AS;   write(nur,'(a)') TRIM(AS)
       if(AS(1:3).eq.'CSF') Exit
      End do

      Call Def_jblocks

      Do jb=1,njbl; ic1=JTc1(jb); ic2=JTc2(jb); nc=ic2-ic1+1

       write(*,*) ic1,ic2,nc,ncfg

       C(1:nc) = -WC(ic1:ic2)
       Call SORTR(nc,C,IP)

       Do j=1,nc; i=IP(j); ic=i+ic1-1
        if(WC(ic).lt.eps_c) Cycle
        Call Print_conf_jj (nur,ic,WC(ic)*cc)
       End do

       write(nur,'(a)') ' *'

      End do

      Backspace(nur)
      write(nur,'(a)') '* '


      End !  Utility jfile


