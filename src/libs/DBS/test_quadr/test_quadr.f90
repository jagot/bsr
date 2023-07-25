!======================================================================
      PROGRAM test_quadrd
!======================================================================
!     Test of one-electron integrals
!----------------------------------------------------------------------
      Use DBS_nuclear, Z => atomic_number, A => atomic_weight
      Use DBS_grid
      Use DBS_galerkin

      Implicit real(8) (A-H,O-Z)

      z = 1.d0;  Call Read_rarg('z',z) 
  
! ... define B-splines:

      Call def_grid(' ',' ',z,1.d0)
      Call alloc_DBS_gauss

      ipri=1; Open(ipri,file='test_quadr.log')

! ... defines the set of hydrogenic orbitals for given Z and all n <= nn
     
      nn = 5;    Call Read_iarg('nn',nn) 

      Call define_orbitals(nn)
    
      Call test_quadr(z)

      END PROGRAM test_quadrd


!========================================================================
      SUBROUTINE define_orbitals(nn)
!========================================================================
! ... define the set of hydrogenic orbitals for n <= nn
!------------------------------------------------------------------------
      Use DBS_nuclear, Z => atomic_number, A => atomic_weight
      Use DBS_grid
      USE DBS_orbitals_pq

      IMPLICIT REAL(8) (A-H,O-Z)

      Do n=1,nn
       Do k = -n,n-1
        if(k.eq.0) Cycle
        i = Ifind_bsorb(n,k,0,2)
        CALL bdcwf_pq(n,k,z,pq(1,1,i),pq(1,2,i)); mbs(i)=ns
       End do
      End do

      nub=1; open(nub,file='dcwf.bsw',form='UNFORMATTED')
      Do i=1,nbf
       write(nub) grid_type,ns,ks
       write(nub) ebs(i),mbs(i)
       write(nub) pq(1:mbs(i),1,i)
       write(nub) pq(1:mbs(i),2,i)
      End do
      Close(nub)

      END SUBROUTINE define_orbitals



!========================================================================
      SUBROUTINE test_quadr(z)
!========================================================================
!     tests the accuracy of one-electron radial integrals on the basis
!     of hydrogenic orbitals. Results are output in file 'quadr_out'.
!------------------------------------------------------------------------
      USE DBS_orbitals_pq

      IMPLICIT NONE
      REAL(8), INTENT(in) :: Z
      REAL(8), EXTERNAL :: QUADR_PQ, QUADR
      INTEGER, EXTERNAL :: Ifind_bsorb
      INTEGER :: iout, i, j, icase
      REAL(8) :: a,ar1,ar2,ar3,am1,am2,am3,br1,br2,br3,bm1,bm2,bm3

      iout=2; Open(iout,file='quadr_out')

!----------------------------------------------------------------------
!                                            accuracy of orthogonality:

      write(iout,'(/a/)')   '  accuracy of orthogonality'
      Do i=1,nbf
       Do j=i,nbf

       if(kbs(i).ne.kbs(j)) Cycle

       a = QUADR_PQ (i,j,0); if(i.eq.j) a = a - 1.d0

       write(iout,'(3x,a,3x,a,D12.3)') EBS(i),EBS(j),a

       End do
      End do

!-----------------------------------------------------------------------
!                                           accuracy of avarages <r^m>:

      write(iout,'(/a/)')   '  accuracy of <r^m>'
      write(iout,'(/4x,5(7x,a5)/)')  &
                 ' <r> ','<r^2>','<r-1>','<r-2>','<r-3>'

      Do i=1,nbf
       write(iout,*)
       Call DCME(nbs(i),kbs(i),z,ar1,ar2,am1,am2,am3)
       br1 = QUADR_PQ(i,i, 1)-ar1
       br2 = QUADR_PQ(i,i, 2)-ar2
       bm1 = QUADR_PQ(i,i,-1)-am1
       bm2 = QUADR_PQ(i,i,-2)-am2
       bm3 = 0.0; if(kbs(i).ne.-1) bm3 = QUADR_PQ(i,i,-3)-am3
       if(ar1.ne.0.d0) br1=br1/ar1
       if(ar2.ne.0.d0) br2=br2/ar2
       if(am1.ne.0.d0) bm1=bm1/am1
       if(am2.ne.0.d0) bm2=bm2/am2
       if(am3.ne.0.d0) bm3=bm3/am3
       write(iout,'(a5,1x,5D12.3)') ebs(i),ar1,ar2,am1,am2,am3
       write(iout,'(6x,   5D12.3)')        br1,br2,bm1,bm2,bm3
      end do

      write(iout,'(/a/)')   '  accuracy of <r^m>'
      write(iout,'(/4x,5(7x,a5)/)')  &
                 ' <r> ','<r^2>','<r-1>','<r-2>','<r-3>'

      Do i=1,nbf
       write(iout,*)
       Call DCME(nbs(i),kbs(i),z,ar1,ar2,am1,am2,am3)
       br1 = QUADR(pq(1,1,i),pq(1,1,i),1)  !-ar1
       br2 = QUADR(pq(1,1,i),pq(1,1,i), 2) !-ar2
       bm1 = QUADR(pq(1,1,i),pq(1,1,i),-1) !-am1
       bm2 = QUADR(pq(1,1,i),pq(1,1,i),-2) !-am2
       bm3 = 0.0; if(kbs(i).ne.-1) bm3 = QUADR(pq(1,1,i),pq(1,1,i),-3) !-am3
       if(ar1.ne.0.d0) br1=br1/ar1
       if(ar2.ne.0.d0) br2=br2/ar2
       if(am1.ne.0.d0) bm1=bm1/am1
       if(am2.ne.0.d0) bm2=bm2/am2
       if(am3.ne.0.d0) bm3=bm3/am3
       write(iout,'(a5,1x,5D12.3)') ebs(i),ar1,ar2,am1,am2,am3
       write(iout,'(6x,   5D12.3)')        br1,br2,bm1,bm2,bm3
      end do

      END SUBROUTINE test_quadr

