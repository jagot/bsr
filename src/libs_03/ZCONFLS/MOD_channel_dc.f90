!===========================================================================
      Module channel_dc
!===========================================================================
!     define scattering channels for one partial wave
!     for the double-continuum case
!---------------------------------------------------------------------------
      Implicit none

      Integer :: lpar =  0 !  total L for given patial wave
      Integer :: ispar=  0 !  (2*S+1) for given patial wave
      Integer :: ipar =  0 !  parity 
      Integer :: nch  =  0 !  number of channels
      Integer :: nch1 =  0 !  number of one-electron channels
      Integer :: nch2 =  0 !  number of two-electron channels
      Integer :: mch  =  0 !  max. number of channels
      Integer :: imch = 64 !  initial prediction of mch

      Integer, allocatable :: iptar(:),ipconf(:),lch1(:),lch2(:)
      Integer, allocatable :: chsym(:),chL(:),chS(:)
      CHARACTER(4), allocatable :: ELC1(:),ELC2(:)
 	
      !   iptar -  pointer on the target state
      !   lch   -  small l for given channel
      !   ELC   -  spectroscopic symbol for given channel
      !   ipconf-  pointer on the last configuration for this channel
	
      Integer :: ncp	     !  number of configurations in perturber
      Integer :: nwp        !  number of orbitals in perturber
      CHARACTER(20) :: AFP,BFP !  file-name for perturber

      End Module channel_dc

    
!=======================================================================    
      Subroutine Alloc_channel_dc(m)
!=======================================================================    
!     allocate arrays in the module "channel_dc"
!-----------------------------------------------------------------------    
      Use channel_dc

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: ia(:)
      Character(4), allocatable :: aa(:)

      if(m.le.0) then
       if(Allocated(iptar)) &
        Deallocate(iptar,ipconf,lch1,lch2,ELC1,ELC2,chsym,chL,chS)
       mch = 0; nch = 0
      elseif(.not.allocated(iptar)) then
       mch = m; nch = 0
       Allocate(iptar(mch),ipconf(0:mch),lch1(mch),lch2(mch), &
        ELC1(mch),ELC2(mch),chsym(mch),chL(mch),chS(mch))
       ipconf=0; lch1=-1; lch2=-1; ELC1 = '    '; ELC2 = '    '
       chL=-1; chS=-1
      elseif(m.le.mch) then
       Return
      elseif(nch.eq.0) then
        Deallocate(iptar,ipconf,lch1,lch2,ELC1,ELC2,chsym,chL,chS)
       mch = m
       Allocate(iptar(mch),ipconf(0:mch),lch1(mch),lch2(mch), &
        ELC1(mch),ELC2(mch),chsym(mch),chL(mch),chS(mch))
       ipconf=0; lch1=-1; lch2=-1; ELC1 = '    '; ELC2 = '    '
       chL=-1; chS=-1
      else
       mch = m
       Allocate(ia(nch))
       ia=iptar(1:nch); Deallocate(iptar)
       Allocate(iptar(mch)); iptar(1:nch)=ia
       ia=ipconf(1:nch); Deallocate(ipconf)
       Allocate(ipconf(mch)); ipconf(1:nch)=ia
       ia=lch1(1:nch); Deallocate(lch1)
       Allocate(lch1(mch)); lch1(1:nch)=ia
       ia=lch2(1:nch); Deallocate(lch2)
       Allocate(lch2(mch)); lch2(1:nch)=ia
       ia=chsym(1:nch); Deallocate(chsym)
       Allocate(chsym(mch)); chsym(1:nch)=ia
       ia=chL(1:nch); Deallocate(chL)
       Allocate(chL(mch)); chL(1:nch)=ia
       ia=chS(1:nch); Deallocate(chS)
       Allocate(chS(mch)); chS(1:nch)=ia
       Deallocate(ia)
       Allocate(aa(nch))
       aa=ELC1(1:nch); Deallocate(ELC1)
       Allocate(ELC1(mch)); ELC1(1:nch)=aa
       aa=ELC2(1:nch); Deallocate(ELC2)
       Allocate(ELC2(mch)); ELC2(1:nch)=aa
       Deallocate(aa)
      end if   

      End subroutine Alloc_channel_dc


!======================================================================
      Subroutine Read_channel_dc(nut,klsp)
!======================================================================
!     read from file 'nut' information for channel klsp 
!----------------------------------------------------------------------
      Use channel_dc
      
      Implicit none
      Integer, intent(in) :: nut,klsp
      Character(20) :: AF
      Character(80) :: line
      Integer :: i,j,n,k,m, nlsp
 
      nlsp = 0
      Call Read_ipar(nut,'nlsp',nlsp); read(nut,*)
      if(nlsp.le.0) Stop 'Read_channel_dc: nlsp <= 0 '
      if(klsp.gt.nlsp) Stop 'R_channel_dc: klsp > nlsp '
      
      Do i = 1,klsp
       read(nut,*) AF,lpar,ispar,ipar,BFP,AFP,ncp,nwp 
      End do

      Call Read_apar(nut,'channels',AF); read(nut,*)

      nch=0
      Do
       read(nut,'(a)') line; if(line(1:3).ne.'par') Cycle
       read(line(14:16),*) i;  if(i.ne.klsp) Cycle
       read(line(31:),*) nch1,nch2
       Call Alloc_channel_dc(nch1+nch2); nch=nch1+nch2
       Do i = 1,nch
        read(nut,*) AF,ELC1(i),ELC2(i),iptar(i),ipconf(i),chL(i),chS(i)
        Call EL4_nlk(ELC1(i),n,lch1(i),k)
        Call EL4_nlk(ELC2(i),n,lch2(i),k)
        chsym(i)=0
        if(lch1(i).eq.lch2(i)) then
         m = (ispar-1)/2 + lpar; chsym(i)=(-1)**m
        end if
       End do
       Exit
      End do

      End Subroutine Read_channel_dc

