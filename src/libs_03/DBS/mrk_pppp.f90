!======================================================================
      Subroutine mrk_pppp(k)
!======================================================================
!     Defines matrix of Rk integrals in the B-spline basis
!     by cell algorithm between P-functions
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_moments
      Use DBS_integrals
      Use Timer

      Implicit none
      Integer, intent(in) :: k
      Integer :: met

! ... check the need of calculations

      if(itype == 'pppp' .and. krk == k) Return

      if (associated(rkb)) nullify(rkb)

      met = -1
      if(ntype.gt.0 .and. k.ge.kra_min .and. k.le.kra_max) then
       rkb => rka(:,:,:,:,k,1)
       if(irka(k,1) == 1) then
        krk=k; itype = 'pppp'; Return
       end if
       met = 0
      end if

      if(met.eq.-1) then
       if(ntype1.eq.0) Call alloc_Rk_integral(ns,ks)
       rkb => rka1(:,:,:,:)
       if(itype1.eq.'pppp'.and.krk1.eq.k) then
        krk=k; itype = 'pppp'; Return
       end if
       met = 1
      end if

      call TimerStart('mrk_pppp')

! ... compute the spline moments:

      call TimerStart('mrk_pppp: moments')
      Call moments_pp(  k   ,kk,nv,rkd1)
      Call moments_pp(-(k+1),kk,nv,rkd2)
      call TimerStop('mrk_pppp: moments')
      call TimerStart('mrk_pppp: diag')
      Call diag_pppp(k)
      call TimerStop('mrk_pppp: diag')

! ... generate the rkb array

      call mrk_common_gen_array(ksp, ksp, rkd1, rkd2, rkd1, rkd2)
#ifdef DEBUG_SPEEDUPS
      call mrk_pppp_gen_array_old()
#endif

      if(met.eq.0) irka(k,1)=1
      if(met.eq.1) then; krk1=k; itype1 = 'pppp'; end if
      krk=k; itype = 'pppp'

      call TimerStop('mrk_pppp')

      End Subroutine mrk_pppp


!======================================================================
      Subroutine diag_pppp(k)
!======================================================================
!     Controls the scaling propeties for diagonal B-spline Rk-interals
!     (not implimented yet)
!----------------------------------------------------------------------
      Use DBS_grid, only: nv

      Implicit none
      Integer :: k,iv

      Do iv=1,nv; Call triang_pppp(k,iv); End do

      END Subroutine diag_pppp


!======================================================================
    Subroutine triang_pppp (k,iv)
!======================================================================
!   Returns the two-dimensional array of B-spline integrals
!               <B_i B_j|r^k/r^(k+1)|B_i' B_j'>
!   over the given triangle diagonal cell
!
!   On entry   iv  -  index of the diagonal cell
!   --------
!
!   On exit    rkd(.,.,iv) - arrays of Rk B-spline integrals for given
!   --------                 interval iv in the reduced-dimension mode
!
!   Calls:   gauleg, zbsplvd
!----------------------------------------------------------------------
    Use DBS_grid
    Use DBS_gauss
    Use DBS_moments
    Use Timer

    Implicit none
    Integer, intent(in) :: k,iv
    Integer :: i,j, ip,jp, ii,jj, m, left,  ik
    Real(8) :: xbase
    Real(8) :: x(ks),w(ks),bi(ks)
    Real(8) :: bspTmp(ks,ksp)
    Real(8) :: INTP(ksp,ksp,ks)
    Real(8) :: a(ksp*(ksp+1)/2,ksp*(ksp+1)/2)

    left=iv+ksp-1;  xbase=tp(left)

!.. setup the gaussian points

    Call gauleg(0.d0,1.d0,x,w,ks)

    Call TimerStart("mrk_pppp: First integration")
    DO m=1,ks

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks
       Call zbsplvd(nsp,ksp,nv, tp,left, 1,gx(i),1,dbip)
       bspTmp(i,1:ksp)= dbip(1,1:ksp,1)
      END DO

! .. and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:)
      if(k>1) then;            gx(:) = gw(:)*gx(:)**k
      else if(k==1) then;      gx(:) = gw(:)*gx(:)
      else if(k==0) then;      gx(:) = gw(:)
      end if

!            / r(iv,m)                             k
! .. INT =  |      bsp(iv,:,j)(r) bsp(iv,:,jp)(r) r  dr
!           / r_iv

      Do j=1,ksp
       gw(:) = gx(:)*bspTmp(:,j)
       Do jp=j,ksp; INTP(j,jp,m)= SUM(gw(:)*bspTmp(:,jp));  END DO
      End do

    END DO    !  over m
    Call TimerStop("mrk_pppp: First integration")

! .. second integration
    Call TimerStart("mrk_pppp: Second integration")
    if(k/=0) then;   gx(:) = grw(iv,:)*grm(iv,:)**(k+1)
    else;            gx(:) = grw(iv,:)*grm(iv,:)
    end if

    ii = 0;  DO i=1,ksp;  DO ip=i,ksp;  ii = ii+1

             bi(:) = pbsp(iv,:,i)*pbsp(iv,:,ip)*gx(:)

    jj = 0;  DO j=1,ksp;  DO jp=j,ksp;  jj = jj+1

             a(ii,jj) =  SUM(bi(:)*INTP(j,jp,:))

             END DO; END DO
             END DO; END DO

    ik = ksp*(ksp+1)/2;  rkd(1:ik,1:ik,iv) = a + TRANSPOSE(a)
    Call TimerStop("mrk_pppp: Second integration")

    End Subroutine triang_pppp


#ifdef DEBUG_SPEEDUPS
    subroutine mrk_pppp_gen_array_old()
      Use DBS_grid
      Use DBS_moments
      Use DBS_integrals
      Use Timer

      Implicit none
      Integer :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
      Real(8) :: c

      Real(8), allocatable :: rkb_ref(:,:,:,:)
      Real(8) :: discrepancy
      Real(8), parameter :: tolerance = sqrt(epsilon(1.d0))

      allocate(rkb_ref(1:ns,1:ns,1:ks,1:ks))
      call TimerStart('mrk_pppp: generate array old')
      rkb_ref=0.d0

      DO jv=1,nv;    jj = 0
      DO jh = 1,ksp; j  = jv  + jh - 1
      DO jhp=jh,ksp; jp = jhp - jh + 1
                     jj = jj  + 1

      DO iv=1,nv;    ii = 0
      DO ih=  1,ksp; i  = iv  + ih - 1
      DO ihp=ih,ksp; ip = ihp - ih + 1
                     ii = ii  + 1

          if     ( iv < jv ) then;   c = rkd1(ii,iv)*rkd2(jj,jv)
          else if( iv > jv ) then;   c = rkd1(jj,jv)*rkd2(ii,iv)
          else;                      c = rkd(ii,jj,iv)
          end if

          rkb_ref(i,j,ip,jp) = rkb_ref(i,j,ip,jp) +  c

      END DO;  END DO;  END DO
      END DO;  END DO;  END DO
      call TimerStop('mrk_pppp: generate array old')

      discrepancy = sum(abs(rkb(1:ns,1:ns,1:ks,1:ks) - rkb_ref(1:ns,1:ns,1:ks,1:ks)))
      write(*,'("Discrepancy: ",e26.16,", tolerance: ",e26.16)') discrepancy, tolerance
      if(discrepancy > tolerance) error stop "Discrepancy exceeds tolerance"

      deallocate(rkb_ref)
    end subroutine mrk_pppp_gen_array_old
#endif
