!====================================================================    
      MODULE spline_Rk_integrals
!====================================================================
!
!     contains the B-spline representation of two-electron integral
!     rkb(i,i';j,j') in symmetric lower-band column storage mode:
!
!            rkb(1:ns,1:ks; 1:ns,1:ks) 
!
!     krk   - multipole index for the integral
!
!--------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(4) :: krk = -2
      REAL(8), ALLOCATABLE :: rkb(:,:,:,:)

      INTEGER(4) :: irkb1
      REAL(8), ALLOCATABLE :: rkb1(:,:,:)
    
      INTEGER(4) :: irkbd,jrkbd, irkbe, jrkbe
      REAL(8), ALLOCATABLE :: rkbd(:,:),rkbe(:,:)

      Integer(4) :: irkbv1, irkbv2, irkbv3
      REAL(8), ALLOCATABLE :: rkbv(:)

      END MODULE spline_Rk_integrals


!====================================================================    
      SUBROUTINE allocate_RK_integrals
!====================================================================    
!     allocates space for spline integrals
!--------------------------------------------------------------------

      Use spline_param
      USE spline_RK_integrals

      Allocate (rkb(ns,ks,ns,ks), rkb1(ns,ns,ks), &
                rkbd(ns,ks), rkbe(ns,ns), rkbv(ns))

      krk = -1; irkb1=0; irkbd=0; jrkbd=0; irkbe=0; jrkbe=0 
                irkbv1=0; irkbv2=0; irkbv3=0      

      END SUBROUTINE allocate_RK_integrals


!====================================================================    
      SUBROUTINE dealloc_RK_integrals
!====================================================================    
!     deallocate arrays in module "spline_Rk_integrals"
!--------------------------------------------------------------------

      USE spline_RK_integrals

      if(krk.le.-1) Return 

      Deallocate (rkb, rkb1, rkbd, rkbe, rkbv)

      krk=-2; irkb1=0; irkbd=0; jrkbd=0; irkbe=0; jrkbe=0 
              irkbv1=0; irkbv2=0; irkbv3=0

      END SUBROUTINE dealloc_RK_integrals


!======================================================================
      SUBROUTINE mrk_cell_lb(k)
!======================================================================
!
!     Defines matrix of Rk integrals in the B-spline basis
!     by cell algorithm
!
!     Calls: rk_moments
!
!----------------------------------------------------------------------
!
!     on entry      k        multipole index
!     --------
!       
!     on exit       rkb     four-dimensional array of Slater integrals 
!     -------               of power k in the B-spline basis
!                           (in module spline-integrals)
!----------------------------------------------------------------------

      USE spline_param
      USE spline_Rk_integrals
      USE spline_moments
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: k
   
      ! .. local variables
   
      INTEGER(4) :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
      REAL(8) :: c
   
      ! .. check the need of calculations
   
      if(krk.eq.k) Return
      if(krk.le.-2) Call allocate_Rk_integrals
   
      ! .. compute the moments in the spline basis
   
      CALL rk_moments(k)
   
      ! .. generate rkb array
   
      rkb=0.d0
   
      DO jv=1,nv;        jj=0
       DO jh=1,ks;       
        DO jhp=jh,ks;    jp=jh-jhp+ks; j=jv+jhp-1; jj=jj+1
   
         DO iv=1,nv;     ii=0
          DO ih=1,ks;    
           DO ihp=ih,ks; ip=ih-ihp+ks; i=iv+ihp-1; ii=ii+1
   
            if( iv < jv ) then;        c = rkd1(ii,iv)*rkd2(jj,jv)
            else if( iv > jv ) then;   c = rkd1(jj,jv)*rkd2(ii,iv)
            else;                      c = rkd(ii,jj,iv) 
            end if
         
            rkb(i,ip,j,jp) = rkb(i,ip,j,jp) +  c 
          
         END DO; END DO; END DO
      END DO; END DO; END DO
   
      krk=k; irkb1=0; irkbd=0; jrkbd=0; irkbe=0; jrkbe=0 
             irkbv1=0; irkbv2=0; irkbv3=0

      END SUBROUTINE mrk_cell_lb


!======================================================================
      Subroutine make_rkb1(k,io)
!======================================================================
!
!     convolutes the rkb(i,i',j,j') array of spline integrals
!     with vector pbs(:,i) from module "spline_orbitals" 
!
!     results in array rkb1(1:ns,1:ns,1:ks) 
!
!     Due to symmetry of rkb array, it does not matter what variable
!     is integrated
!
!----------------------------------------------------------------------

      USE spline_param
      USE spline_Rk_integrals
      Use spline_orbitals

      IMPLICIT NONE

      INTEGER(4) :: k,io, j,jth

      if(irkb1.eq.io) Return 

      if(k.ne.krk) Call mrk_cell_lb(k) 

      do jth=1,ks;  do j=ks+1-jth,ns
       Call bxv(ks,ns,rkb(1,1,j,jth),pbs(1,io),rkb1(1,j,jth))
      end do; end do

      irkb1 = io

      End Subroutine make_rkb1
       

!======================================================================
      Subroutine make_rkd(k,io,jo)
!======================================================================
!
!     convolutes the rkb(i,j;i',j') array of spline integrals
!     with matrix ds(i,j) from module "spline_densities" 
!
!     results in array rkbd(1:ns,1:ks) 
!
!     Due to symmetry of rkb array, it does not matter what variable
!     is integrated,  (i,j) or (i',j') 
!
!----------------------------------------------------------------------

      USE spline_param
      USE spline_Rk_integrals
      USE spline_densities

      IMPLICIT NONE

      INTEGER(4) :: k,io,jo, j,jp

      if(io.eq.irkbd.and.jo.eq.jrkbd) Return 

      if(k.ne.krk) Call mrk_cell_lb(k) 

      Call make_density(io,jo,'s')

      rkbd = 0.d0
      do jp = 1,ks
       do j = ks+1-jp,ns
         rkbd(j,jp) = SUM(ds(:,:)*rkb(:,:,j,jp))
       end do
      end do

      irkbd = io; jrkbd = jo

      End Subroutine make_rkd


!======================================================================
      Subroutine make_rke(k,io,jo)
!======================================================================
!
!     convolutes the rkb(i,j;i',j') array of spline integrals
!     with matrix dx(i,j) from module "spline_densities" 
!
!     results in array rkbe(1:ns,1:ns) 
!
!     Due to symmetry of rkb array, it does not matter what variable
!     is integrated,  (i,j') or (i',j) 
!
!----------------------------------------------------------------------

      USE spline_param
      USE spline_Rk_integrals
      USE spline_densities

      IMPLICIT NONE

      INTEGER(4) :: k,io,jo, i,ip,ith, j,jp,jth
      Real(8) :: c

      if(io.eq.irkbe.and.jo.eq.jrkbe) Return 

      if(k.ne.krk) Call mrk_cell_lb(k) 

      Call make_density(io,jo,'x')

      rkbe = 0.d0

      do jth=1,ks;  do j=ks+1-jth,ns; jp=j+jth-ks
      do ith=1,ks;  do i=ks+1-ith,ns; ip=i+ith-ks
       c = rkb(i,ith,j,jth)
                                   rkbe(i ,jp) = rkbe(i ,jp) + c * dx(j ,ip)
       if(ith.ne.ks)               rkbe(ip,jp) = rkbe(ip,jp) + c * dx(j ,i )
       if(jth.ne.ks)               rkbe(i ,j ) = rkbe(i ,j ) + c * dx(jp,ip)
       if(ith.ne.ks.and.jth.ne.ks) rkbe(ip,j ) = rkbe(ip,j ) + c * dx(jp,i ) 
      
      end do; end do
      end do; end do

      irkbe = io; jrkbe = jo

      End Subroutine make_rke


!======================================================================
      Subroutine make_rkbv(k,i1,i2,i3)
!======================================================================
!
!     convolutes the rkb(i,P1;P2,P3) array of spline integrals
!     over three bound orbitals 
!
!     results in array rkbv(1:ns) 
!
!     Due to symmetry of rkb array, it does not matter what variable
!     is integrated,  (i,j;i',j') 
!
!----------------------------------------------------------------------

      USE spline_param
      USE spline_Rk_integrals
      USE spline_orbitals
      USE spline_galerkin

      IMPLICIT NONE

      INTEGER(4) :: k,i1,i2,i3

      if(i1.eq.irkbv1.and.i2.eq.irkbv2.and.i3.eq.irkbv3) Return 

      if(k.ne.krk) Call mrk_cell_lb(k) 

      Call make_rkd(k,i1,i3)

      Call bxv(ks,ns,rkbd,pbs(1,i2),rkbv)

      irkbv1=i1; irkbv2=i2; irkbv3=i3

      End Subroutine make_rkbv
