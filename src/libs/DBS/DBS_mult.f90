!====================================================================
      Module DBS_mult
!====================================================================
!     contains matrices for  the radiative transition integrals
!     in B-spline basis (kp and kq)
!--------------------------------------------------------------------
!     long wavefunction case:
!--------------------------------------------------------------------
!
!     dipL_pp => Int[ {P} * r^k     * {P},  dr]
!     dipL_QQ => Int[ {Q} * r^k     * {Q},  dr]
!     dipV_PQ => Int[ {P} * r^(k-1) * {Q},  dr]
!
!--------------------------------------------------------------------
!     relativistic case:
!--------------------------------------------------------------------
!
!     dipk_pp  =>  Int[ P * j_k (omega*r/c) * P,  dr]
!     dipk_qq  =>  Int[ Q * j_k (omega*r/c) * Q,  dr]
!
!     dipkm_pq =>  Int[ P * j_(k-1)(omega*r/c) * Q,  dr]
!     dipkp_pq =>  Int[ P * j_(k+1)(omega*r/c) * Q,  dr]
!
!--------------------------------------------------------------------
      Implicit none

      Integer :: kmult = -1    ! multipole index for current save
      Real(8) :: wmult = 0.d0  ! omega factor

      Real(8), allocatable :: dipL_pp(:,:),dipL_qq(:,:),dipV_pq(:,:)

! ... relativistic case (wmult <> 0):

      Real(8), allocatable :: dipk_pp(:,:),dipk_qq(:,:)

      Real(8), allocatable :: dipkm_pq(:,:),dipkp_pq(:,:)

      End Module DBS_mult


!======================================================================
      Subroutine dip_mat_pq(wfact,kpol)
!======================================================================
!     prepare B-spline dipole matrixes
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_mult

      Implicit none
      Real(8), intent(in) :: wfact
      Integer, intent(in) :: kpol
      Real(8), allocatable :: yb(:,:,:),bes(:)
      Integer :: i,j

      if(kpol.eq.kmult.and.wfact.eq.wmult) Return

! ... long wavelength case:

      if(wfact.eq.0.d0) then

       if(.not.Allocated(dipL_pp)) Allocate(dipL_pp(ns,ns))
       if(.not.Allocated(dipL_qq)) Allocate(dipL_qq(ns,ns))
       if(.not.Allocated(dipV_pq)) Allocate(dipV_pq(ns,ns))

       Call ZINTYk (kpol,  ksp,ksp,pbsp,pbsp,ns,dipL_pp)
       Call ZINTYk (kpol,  ksq,ksq,qbsp,qbsp,ns,dipL_qq)
       Call ZINTYk (kpol-1,ksp,ksq,pbsp,qbsp,ns,dipV_pq)

      else

! ... relativistic case:

      if(.not.Allocated(dipk_pp )) Allocate(dipk_pp (ns,ns))
      if(.not.Allocated(dipk_qq )) Allocate(dipk_qq (ns,ns))
      if(.not.Allocated(dipkm_pq)) Allocate(dipkm_pq(ns,ns))
      if(.not.Allocated(dipkp_pq)) Allocate(dipkp_pq(ns,ns))
      Allocate(yb(nv,ks,0:kpol+1),bes(0:kpol+1))

      Do i=1,nv; Do j=1,ks
       Call SPBES(wfact*gr(i,j),kpol+2,bes(0))
       yb(i,j,:) = bes(:) * grw(i,j)
      End do; End do

      Call ZINTYm (nv,ks,ksp,ksp,pbsp,pbsp,yb(1,1,kpol  ),ns,dipk_pp )
      Call ZINTYm (nv,ks,ksq,ksq,qbsp,qbsp,yb(1,1,kpol  ),ns,dipk_qq )
      Call ZINTYm (nv,ks,ksp,ksq,pbsp,qbsp,yb(1,1,kpol-1),ns,dipkm_pq)
      Call ZINTYm (nv,ks,ksp,ksq,pbsp,qbsp,yb(1,1,kpol+1),ns,dipkp_pq)

      Deallocate(yb,bes)

      end if

      wmult=wfact; kmult=kpol

      End Subroutine dip_mat_pq


!======================================================================
      Real(8) Function JKDIP_pq(i,j,kpol,wfact)
!======================================================================
!     J_k = Int [ (P_i(r)*P_j(r)+Q_i(r)*Q_j(r)) * j_k(omega) , dr]
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_mult
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: i,j,kpol
      Real(8), intent(in) :: wfact
      Real(8), external :: vav

      Call dip_mat_pq(wfact,kpol)

      JKDIP_pq = vav (nsp,ksp,nsp,ksp,ns,dipk_pp,pq(1,1,i),pq(1,1,j),'x') +  &
                 vav (nsq,ksq,nsq,ksq,ns,dipk_qq,pq(1,2,i),pq(1,2,j),'x')

      End Function JKDIP_pq


!======================================================================
      Subroutine IKDIP_pq(i,j,kpol,wfact,CP,CM,met)
!======================================================================
!     IK+- = Int [ (P_i(r) * j_k(omega)* Q_j(r)), dr] +-
!            Int [ (Q_i(r) * j_k(omega)* P_j(r)), dr]
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_mult
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: i,j,kpol,met
      Real(8), intent(in) :: wfact
      Real(8) :: CP,CM,S1,S2,v(ns)
      Real(8), external :: vav

      Call dip_mat_pq(wfact,kpol)

      Select case(met)
       Case(-1)

        S1 = vav (nsp,ksp,nsq,ksq,ns,dipkm_pq,pq(1,1,i),pq(1,2,j),'x')
        S2 = vav (nsp,ksp,nsq,ksq,ns,dipkm_pq,pq(1,1,j),pq(1,2,i),'x')

       Case(+1)

        S1 = vav (nsp,ksp,nsq,ksq,ns,dipkp_pq,pq(1,1,i),pq(1,2,j),'x')
        S2 = vav (nsp,ksp,nsq,ksq,ns,dipkp_pq,pq(1,1,j),pq(1,2,i),'x')

       Case Default

        Stop 'IKDIP_pq: unknown mode'

      End Select
      CP = S1 + S2
      CM = S1 - S2

      End Subroutine IKDIP_pq


!======================================================================
      Real(8) Function dipM_pq(i,j,kpol,wfact)
!======================================================================
!     D_m = (2k+1) / sqrt[k*(k+1)] * (kappa_i+kappa_j) * I_k^+
!----------------------------------------------------------------------
      Use DBS_mult
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: i,j,kpol
      Real(8), intent(in) :: wfact
      Real(8) :: C,CP,CM
      Integer, external :: ITRA
      Real(8), external :: DQUADRM_pq

      dipM_pq = 0.d0

      if(wfact.eq.0.d0) then; dipM_pq=DQUADRM_pq(i,j,kpol); Return; end if

      if(ITRA(jbs(i),kpol+kpol,jbs(j)).eq.0) Return
      if(mod(lbs(i)+kpol+lbs(j),2).ne.1) Return

      C = kpol*(kpol+1)
      C = (kbs(i)+kbs(j))/SQRT(C)*(kpol+kpol+1)
      Call IKDIP_pq(i,j,kpol,wfact,CP,CM,1)
      dipM_pq = CP * C

      End Function dipM_pq


!======================================================================
      Real(8) Function DQUADRM_pq(i,j,k)
!======================================================================
!     magnetic dipole radial integral in low frequency limit
!
!     (kj+ki)/(k+1) (<Pi|r^k|Qj>+<Qi|r^k|Pj>)
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_mult
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: i,j,k
      Real(8) :: spq, sqp
      Integer, external :: ITRA
      Real(8), external :: vav

      DQUADRM_pq = 0.d0
      if(k.lt.1) Stop 'DQUADRM: kpol < 1'

      if(ITRA(jbs(i),k+k,jbs(j)).eq.0) Return
      if(mod(lbs(i)+k+lbs(j),2).ne.1) Return

      spq = vav (nsp,ksp,nsq,ksq,ns,dipkp_pq,pq(1,1,i),pq(1,2,j),'x')
      sqp = vav (nsp,ksp,nsq,ksq,ns,dipkp_pq,pq(1,1,j),pq(1,2,i),'x')

      DQUADRM_pq = (spq + sqp)*(kbs(i)+kbs(j))/(k+1)

      End function DQUADRM_pq


!======================================================================
      Subroutine dipE_pq(i,j,kpol,wfact,CL,CV)
!======================================================================
!     Define:  kk = kappa_i-kappa_j
!              II+ = kk * I+_(k+1) + (k+1) * I-_(k+1)
!              II- = kk * I+_(k-1) -  k    * I-_(k-1)
!
! Coulomb (length) gauge:  D = 1/sqrt(k)/sqrt(k+1) [k*II+  -  (k+1)*II-]
! Babushkin (velocity)  :  D = (2k+1)/sqrt(k)/sqrt(k+1) [II+  +- Jk]
!-----------------------------------------------------------------------
      Use DBS_orbitals_pq, only: kbs,lbs,jbs,ebs

      Implicit none
      Integer, intent(in) :: i,j,kpol
      Real(8), intent(in) :: wfact
      Real(8) :: C,CP1,CM1,CP2,CM2,CJ,CL,CV,C1,C2
      Integer :: k,li,lj,ji,jj
      Integer, external :: ITRA
      Real(8), external :: DQUADRL_pq,DQUADRV_pq,JKDIP_pq

      CL = 0.d0; CV = 0.d0
      if(ITRA(jbs(i),kpol+kpol,jbs(j)).eq.0) Return
      if(mod(lbs(i)+kpol+lbs(j),2).ne.0) Return

      if(wfact.eq.0.d0) then
       CL = DQUADRL_pq(i,j,kpol)
       CV = DQUADRV_pq(i,j,kpol)
       Return
      end if

      Call IKDIP_pq(i,j,kpol,wfact,CP1,CM1,-1)
      Call IKDIP_pq(i,j,kpol,wfact,CP2,CM2, 1)
      CJ = JKDIP_pq(i,j,kpol,wfact) * (kpol+kpol+1)

      C = kpol+1; C = sqrt(C/kpol)
      k = kbs(i)-kbs(j)
      C1 = (k*CP2+(kpol+1)*CM2)
      C2 = (k*CP1- kpol   *CM1)
      CV = C1/C - C*C2
      CL = CV + C*(C1+C2-CJ)              !  ??? + or -

      End Subroutine dipE_pq


!======================================================================
      Real(8) Function DQUADRL_pq (i,j,k)
!======================================================================
!     electriv dipole radial integral in low frequency limit
!     in  Babushkin (length) gauge:   <Pi|r^k|Pj> + <Qi|r^k|Qj>
!----------------------------------------------------------------------
      Use DBS_orbitals_pq, only: kbs,lbs,jbs

      Implicit none
      Integer, intent(in) :: i,j,k
      Integer :: li,lj,ji,jj
      Integer, external :: l_kappa,j_kappa,ITRA
      Real(8), external :: QUADR_pq

      DQUADRL_pq = 0.d0
      if(k.lt.1) Stop 'DQUADRL_pq: kpol < 1'

      if(ITRA(jbs(i),k+k,jbs(j)).eq.0) Return
      if(mod(lbs(i)+k+lbs(j),2).ne.0) Return

      DQUADRL_pq = QUADR_pq(i,j,k)

      End Function DQUADRL_pq


!======================================================================
      Real(8) Function DQUADRV_pq(i,j,k)
!======================================================================
!     electriv dipole radial integral in low frequency limit
!     in Coulomb (velocity) gauge:
!       (ki-kj-k)*<Pi|r^k-1|Qj> + (ki-kj+k)*<Qi|r^k-1|Pj>)
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: i,j,k
      Real(8) :: spq, sqp
      Integer, external :: ITRA
      Real(8), external :: quadrm

      DQUADRV_pq = 0.d0
      if(k.lt.1) Stop 'DQUADRV_pq: kpol < 1'

      if(ITRA(jbs(i),k+k,jbs(j)).eq.0) Return
      if(mod(lbs(i)+k+lbs(j),2).ne.0) Return

      spq =  quadrm(nv,ks,ksp,ksq,k-1,pbsp,qbsp,pq(1,1,i),pq(1,2,j))
      sqp =  quadrm(nv,ks,ksq,ksp,k-1,qbsp,pbsp,pq(1,2,i),pq(1,1,j))

      DQUADRV_pq = (kbs(i)-kbs(j)-k)*spq  +  (kbs(i)-kbs(j)+k)*sqp

      End function DQUADRV_pq


