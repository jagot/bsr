!======================================================================
      Module target
!======================================================================
!     contains description of target states
!----------------------------------------------------------------------
      Implicit none

      Integer :: ntarg	 !  number of target states
      Integer :: nphys	 !  number of physical states
      Integer :: nelc	 !  number of atomic electrons
      Integer :: nz      !  atomic number
      Integer :: nct	 !  total number of target config.s 
      Integer :: nwt     !  total number of taget orbitals

      Integer, allocatable :: istarg(:) !  (2*S+1) for target i 
      Integer, allocatable :: ltarg(:)  !  total L 
      Integer, allocatable :: iptarg(:) !  parity (+-1)
      Integer, allocatable :: jtarg(:)  !  (2*J+1) 
      Real(8), allocatable :: etarg(:)  !  target energy in au

      Integer, allocatable :: nctarg(:) !  number of target configurations
      Integer, allocatable :: nwtarg(:) !  number of new orbitals
      Integer, allocatable :: ictarg(:) !  pointer to target i in conf.list

      Character(20), allocatable :: AFT(:) ! given file-names for target states
      Character(20), allocatable :: BFT(:) ! original file-names for target states

      Character(2) :: COUPLING = 'LS' ! LS, JK or JJ coupling mode

      End Module target


!======================================================================
      Subroutine allocate_target(m)
!======================================================================
!     allocate (deallocate) space in module "target" 
!----------------------------------------------------------------------
      Use target
      Implicit none
      Integer, Intent(in) :: m      

      if(m.le.0) then
       if(allocated(istarg)) Deallocate(istarg,ltarg,iptarg,jtarg, &
                nctarg, nwtarg, ictarg, etarg, AFT, BFT)
       ntarg = 0
      else
       if(allocated(istarg)) Deallocate(istarg,ltarg,iptarg,jtarg, &
                nctarg, nwtarg, ictarg, etarg, AFT, BFT)
       ntarg = m
       Allocate(istarg(ntarg),  ltarg(ntarg), iptarg(ntarg), &
                nctarg(ntarg), nwtarg(ntarg), ictarg(ntarg), &
                etarg(ntarg), AFT(ntarg),BFT(ntarg), jtarg(ntarg) )
      end if

      nct = 0
      nwt = 0

      End Subroutine allocate_target


!======================================================================
      Subroutine R_target(nut)
!======================================================================
!     read from file 'nut' target information 
!----------------------------------------------------------------------
      Use target 
      Implicit none
      Integer, Intent(in) :: nut
      Integer :: i
      Character(20) :: AF
 
      Call Read_ipar(nut,'nelc',nelc)
      Call Read_ipar(nut,'nz',nz)
      Call Read_ipar(nut,'nphys',nphys)
      Call Read_apar(nut,'coupling',coupling)
      Call Read_ipar(nut,'ntarg',ntarg)

      if(ntarg.le.0) Stop 'R_targ: ntarg <= 0 '
      i = ntarg; Call Allocate_target(i)
        
      Call Read_ipar(nut,'nwt',nwt)
      Call Read_ipar(nut,'ntarg',ntarg)

      nct = 0
      read(nut,*) 
      Do i=1,ntarg
       read(nut,*) BFT(i),AFT(i),ltarg(i),istarg(i),iptarg(i),etarg(i), &
                   nctarg(i),nwtarg(i)
       nct=nct+nctarg(i)
!       nwt=nwt+nwtarg(i)
       ictarg(i) = nct
       jtarg(i)=0; if(istarg(i).eq.0) jtarg(i)=ltarg(i)+1 
      End do

      End Subroutine R_target


!======================================================================
      Subroutine Write_target_LS(nut)
!======================================================================
!     write target information to file 'nut'  
!----------------------------------------------------------------------
      Use target
      Implicit none
      Integer, Intent(in) :: nut
      Integer :: i

      rewind(nut)
      write(nut,'(a)') 'TITLE' 
      write(nut,'(80(''-''))')
      write(nut,'(a,i4,5x,a)') &
                'nz    = ',nz,   ' !   nuclear charge' 
      write(nut,'(a,i4,5x,a)') &
                'nelc  = ',nelc, ' !   number of electrons'
      write(nut,'(a,a2,4x,a)') &
                'coupling = ',coupling, ' !  coupling scheme'

      write(nut,'(80(''-''))')
      write(nut,'(a,i4,5x,a)') &
                'ntarg = ',ntarg,' !   number of target states'
      write(nut,'(80(''-''))')
      Do i=1,ntarg
       write(nut,'(a20,2x,a10,2x,3i4,f18.8,2i5)') AFT(i),BFT(i), &
        ltarg(i),istarg(i),iptarg(i),etarg(i),nctarg(i),nwtarg(i)
      End do
      write(nut,'(80(''-''))')
      write(nut,'(a,i7)') 'nct =',nct
      write(nut,'(a,i7)') 'nwt =',nwt
      write(nut,'(80(''-''))')

      End Subroutine Write_target_LS


!======================================================================
      Subroutine Print_target
!======================================================================
! ... print example of target file
!----------------------------------------------------------------------
      Character :: AS(80)

      Call Find_free_unit(nut)
      Open(nut,file='target')

      write(nut,'(a)') & 
'title of the case:  e + ... ',                                            &
'------------------------------------------------------------------------',&
'coupling = LS    !   coupling scheme',                                    &        
'nz =  14         !   nuclear charge',                                     &        
'nelc = 13        !   number of electrons',                                &       
'------------------------------------------------------------------------',&
'ntarg =  7       !   number of target states',                            &       
'------------------------------------------------------------------------',&
'3s2_3p       ',                                                           &
'3s3p2_4P     ',                                                           &
'3s3p2_2D     ',                                                           &
'4s_2S        ',                                                           & 
'3s3p2_2S     ',                                                           &
'3d_2D        ',                                                           &
'3s3p2_2P     ',                                                           & 
'------------------------------------------------------------------------',&
'nlsp = 8         !   number of partial waves',                            &     
'------------------------------------------------------------------------',&        
'001    0    1    1   no       ',                                          & 
'002    0    3    1   no       ',                                          &
'003    1    1    1   no       ',                                          &
'004    1    3    1   3p2_3P   ',                                          &
'005    1    1   -1   3p3_1P   ',                                          &
'006    1    3   -1            ',                                          &
'007    2    1    1            ',                                          &
'008    2    3    1            ',                                          &
'------------------------------------------------------------------------',&        
'kpert = 2       !  number of additional perturbers (optional)',           &                        
'------------------------------------------------------------------------',&        
'001  3p2_1S  ',                                                           &
'006  3p3_3P  ',                                                           &
'------------------------------------------------------------------------'

      End Subroutine Print_target
