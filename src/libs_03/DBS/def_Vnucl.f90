!======================================================================
      Subroutine Def_Vnucl
!======================================================================
!     define nuclear potential for p- and q-splines
!---------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_nuclear

      Implicit none
      Integer :: i,j
      Real(8), external :: V_nuclear

      if(allocated(VR_nucl)) Deallocate(VR_nucl)
      Allocate(VR_nucl(nv,ks))

      Do i=1,nv
       Do j=1,ks
        VR_nucl(i,j) = V_nuclear(gr(i,j)) * grw(i,j)
       End do
      End do

      Call ZINTYm (nv,ks,ksp,ksp,pbsp,pbsp,VR_nucl,ns,fpb_nucl)
      Call ZINTYm (nv,ks,ksq,ksq,qbsp,qbsp,VR_nucl,ns,fqb_nucl)

      End Subroutine Def_Vnucl

!======================================================================
      Subroutine Def_ro_fermi
!======================================================================
!     define nuclear potential for p- and q-splines
!---------------------------------------------------------------------
      Use DBS_nuclear, Z => atomic_number, a => a_fermi, c => c_fermi
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer :: i,j
      Real(8) :: S

      Do i=1,nv
       Do j=1,ks
        ygw(i,j) = one/(one+exp((gr(i,j)-c)/a)) * grw(i,j) * gr(i,j) * gr(i,j)
       End do
      End do
      S =  SUM(ygw(:,:))

      ro_fermi = Z / (S * 4 * pi)

      End Subroutine Def_ro_fermi


!======================================================================
      Subroutine Def_Znucl
!======================================================================
!     define nuclear charge destribution in gausian points
!---------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_nuclear

      Implicit none
      Integer :: i,j
      Real(8), external :: Z_nuclear

      if(allocated(ZR_nucl)) Deallocate(ZR_nucl)
      Allocate(ZR_nucl(nv,ks))

      Do i=1,nv
       Do j=1,ks
        ZR_nucl(i,j) = Z_nuclear(gr(i,j)) * grw(i,j)
       End do
      End do

      End Subroutine Def_Znucl

