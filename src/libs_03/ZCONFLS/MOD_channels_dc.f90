!===========================================================================
      Module channels_dc
!===========================================================================
!     define scattering channels for one partial wave
!---------------------------------------------------------------------------
      Implicit none

      Integer :: nlsp = 0  !  number of partial waves                 
      Integer :: mlsp = 0  !  max.dimension                 

      Integer, allocatable :: lpar (:)    !  total L                      
      Integer, allocatable :: ispar(:)    !  (2*S+1)              
      Integer, allocatable :: ipar (:)    !  parity (+1 | -1)             
      Integer, allocatable :: nch  (:)    !  number of channels           
      Integer, allocatable :: nch1 (:)    !  one-continuum                
      Integer, allocatable :: nch2 (:)    !  two-continuum                
    
      CHARACTER(3), allocatable :: Tpar(:)   !  spectroscopic notation 

      Integer, allocatable :: iptar(:),ipconf(:),lch1(:),lch2(:)
      Integer, allocatable :: chsym(:),chL(:),chS(:)
      CHARACTER(4), allocatable :: ELC1(:),ELC2(:)
 	
      !   iptar -  pointer on the target state
      !   lch   -  small l for given channel
      !   ELC   -  spectroscopic symbol for given channel
      !   ipconf-  pointer on the last configuration for this channel
	
      Integer, allocatable :: ipch(:)     !  partial wave pointer

      !   perturber information

      Integer, allocatable :: ncp(:)
      Integer, allocatable :: nwp(:)
      CHARACTER(20), allocatable :: AFP(:),BFP(:) 
      
      Integer :: mch = 0  !  max.number of channels                 

!      Character(3) :: ALSP,BLSP

      End Module channels_dc


!=======================================================================    
      Subroutine Alloc_channels_dc
!=======================================================================    
!     allocate arrays in the module "channels_dc"
!-----------------------------------------------------------------------        
      Use channels_dc

      if(nlsp.le.0) Stop ' Alloc_channels_dc: nlsp <= 0 '
      if(mlsp.le.0) Stop ' Alloc_channels_dc: mlsp <= 0 '

      Allocate(lpar(nlsp),ispar(nlsp),ipar(nlsp),Tpar(nlsp),&
             nch(nlsp),nch1(nlsp),nch2(nlsp), &
             ipch(nlsp),iptar(mlsp),ipconf(mlsp), &
             chsym(mlsp),chL(mlsp),chS(mlsp), &
             ELC1(mlsp), ELC2(mlsp), lch1(mlsp),lch2(mlsp), &
             AFP(nlsp),BFP(nlsp),ncp(nlsp),nwp(nlsp))

      End subroutine Alloc_channels_dc


!======================================================================
      Subroutine Read_channels_dc(nut)
!======================================================================
!     read from file 'nut' channels information 
!----------------------------------------------------------------------
      Use channels_dc
      
      Integer, Intent(in) :: nut
      Character(20) :: AF
      Character(80) :: line
 
! ... define mlsp:

      nlsp = 0
      Call Read_ipar(nut,'nlsp',nlsp); read(nut,*)

      Call Read_apar(nut,'channels',AF)
      mlsp = 0
      Do i = 1,nlsp
       read(nut,*); read(nut,'(a80)') line
       read(line(30:),*) i1,i2; mlsp = mlsp + i1 + i2
       Do j = 1,i1+i2; read(nut,*);  End do
      End do

! ... define nlsp:

      nlsp = 0
      Call Read_ipar(nut,'nlsp',nlsp); read(nut,*)

      Call Alloc_channels_dc      

      Do i = 1,nlsp
       read(nut,*) Tpar(i),lpar(i),ispar(i),ipar(i), &
                   BFP(i),AFP(i),ncp(i),nwp(i)
      End do

      Call Read_apar(nut,'channels',AF)
 
      ip = 0
      Do i = 1,nlsp
       read(nut,*); read(nut,'(a80)') line
       read(line(31:),*) nch1(i),nch2(i); nch(i) = nch1(i)+nch2(i)
       Do j = 1,nch(i); jp=ip+j
        read(nut,'(a8,a4,2x,a4,4i8)') &
         AF,ELC1(jp),ELC2(jp),iptar(jp),ipconf(jp),chL(jp),chS(jp)
        Call EL4_nlk(ELC1(jp),n,lch1(jp),k)
        lch2(jp) = -1
        if(ELC2(jp).ne.'    ') Call EL4_nlk(ELC2(jp),n,lch2(jp),k)
        chsym(jp)=0
        if(lch1(jp).eq.lch2(jp)) chsym(jp)=(-1)**(chL(jp)+(chS(jp)-1)/2) 
       End do
       ipch(i) = ip; ip = ip + nch(i)
      End do

      mch = 0
      Call Read_ipar(nut,'max_ch',mch)

      End Subroutine Read_channels_dc
