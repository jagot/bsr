!===========================================================================
      Module dcs_pseudo
!===========================================================================
!     contain arrays required for differential cross sections
!     in case of many pseudo-states
!---------------------------------------------------------------------------
      Implicit none

      Integer :: nps  = 0                 !  number of psedo_states
      Integer :: ng   = 0                 !  number of angles
      Integer :: i16  = 16                !  dcs unit
      Integer :: Lmax = 0                 !
      Integer :: Smax = 0                 !

      Real(8) :: EK  = 0.d0               !  initial energy in Ry
      Real(8) :: EA  = 0.d0               !  absolute energy (a.u.)
      Real(8) :: EK1 = 0.d0               !  primery electron energy
      Real(8) :: EK2 = 0.d0               !  secondary electron energy

      Real(8), allocatable :: g(:)        !  angles
      Real(8), allocatable :: dcs(:,:)    !  dcs         (ig,is)
      Real(8), allocatable :: f(:,:,:,:)  !  amplitudes  (k,m,is,ig)
      Real(8), allocatable :: EP(:)       !  energy of state

      Integer, allocatable :: ILPS(:),ISPS(:),IPPS(:)

      Real(8)  :: eps_ek = 1.d-7

      End Module dcs_pseudo


!======================================================================
      Subroutine Read_dcs_pseudo(AF)
!======================================================================
!     read from file AF dcs information for pseudo-states
!----------------------------------------------------------------------
      Use dcs_pseudo

      Implicit none
      Character(*) :: AF
      Character(200) :: AS
      Integer :: nu,i,k,l,m, ig,is
      Real(8) :: s
      Integer, external :: Ifind_position

      Call Check_file(AF)
      Call Find_free_unit(nu)
      open(nu,file=AF)

      EK = 0.d0
      Call Read_rpar(nu,'EK',EK)
      if(EK.eq.0.d0) Stop 'Read_dcs_pseudo: EK ?'

      ng = 0
      Call Read_ipar(nu,'ng',ng)
      if(ng.eq.0) Stop 'Read_dcs_pseudo: ng ?'
      if(allocated(G)) Deallocate(G)
      Allocate(G(ng))
      Do ig=1,ng; read(nu,*) g(ig); End do

      nps = 0
      Call Read_ipar(nu,'nstates',nps)
      if(nps.eq.0) Stop 'Read_dcs_pseudo: nps ?'

      lmax = -1
      Call Read_ipar(nu,'lmax',lmax)
      if(lmax.eq.-1) Stop 'Read_dcs_pseudo: lmax ?'

      if(allocated(DCS)) deallocate(DCS,f,EP,ILPS,ISPS,IPPS)
      Allocate(DCS(ng,nps),F(2,-lmax:lmax,nps,ng), &
               EP(nps),ILPS(nps),ISPS(nps),IPPS(nps))

      rewind(nu)
      i=Ifind_position(nu,'is,E,ilsp')
      Do is = 1,nps
       read(nu,'(a)') AS
       i = INDEX(AS,'=')+1
       read(AS(i:),*) i,EP(is),ILPS(is),ISPS(is),IPPS(is)
       l = ILPS(is)
       Do ig = 1,ng
        read(nu,*) s,dcs(ig,is),((f(k,m,is,ig),k=1,2),m=-l,l)
       End do
      End do
      Close(nu)
      Smax = maxval(ISPS)

      End Subroutine Read_dcs_pseudo

