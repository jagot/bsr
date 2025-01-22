!======================================================================
      Subroutine mrk_pqpq(k)
!======================================================================
!     Defines matrix of Rk integrals in the B-spline basis
!     by cell algorithm for mixed PQ case
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_moments
      Use DBS_integrals
      Use Timer

      Implicit none
      Integer, intent(in) :: k
      Integer :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp, met
      Integer, external :: Icheck_rka
      Real(8) :: c

! ... check the need of calculations

      if(itype == 'pqpq' .and. krk == k) Return

      if (associated(rkb)) nullify(rkb)

      met = -1
      if(ntype.gt.0 .and. k.ge.kra_min .and. k.le.kra_max) then
       rkb => rka(:,:,:,:,k,3)
       if(irka(k,3) == 1) then
        krk=k; itype = 'pqpq'; Return
       end if
       met = 0
      end if

      if(met.eq.-1) then
       if(ntype1.eq.0) Call alloc_Rk_integral(ns,ks)
       rkb => rka1(:,:,:,:)
       if(itype1.eq.'pqpq'.and.krk1.eq.k) then
        krk=k; itype = 'pqpq'; Return
       end if
       met = 1
      end if

      call TimerStart('mrk_pqpq')

! ... compute the spline moments:

      call TimerStart('mrk_pqpq: moments')
      Call moments_pp(  k   ,kk,nv,rkd1)
      Call moments_qq(-(k+1),kk,nv,rkd2)
      Call moments_qq(  k   ,kk,nv,rkd3)
      Call moments_pp(-(k+1),kk,nv,rkd4)
      call TimerStop('mrk_pqpq: moments')
      call TimerStart('mrk_pqpq: diag')
      Call diag_pqpq(k)
      call TimerStop('mrk_pqpq: diag')

! ... generate the rkb array

      call TimerStart('mrk_pqpq: generate array')
      rkb=0.d0

      DO jv = 1,nv;   jj = 0
      DO jh = 1,ksq;  j  = jv  + jh - 1
      DO jhp=jh,ksq;  jp = jhp - jh + 1
                      jj = jj  + 1

      DO iv = 1,nv;   ii = 0
      DO ih = 1,ksp;  i  = iv  + ih - 1
      DO ihp=ih,ksp;  ip = ihp - ih + 1
                      ii = ii  + 1

          if     ( iv < jv ) then;   c = rkd1(ii,iv)*rkd2(jj,jv)
          else if( iv > jv ) then;   c = rkd3(jj,jv)*rkd4(ii,iv)
          else;                      c = rkd(ii,jj,iv)
          end if

          rkb(i,j,ip,jp) = rkb(i,j,ip,jp) +  c

      END DO;  END DO;  END DO
      END DO;  END DO;  END DO

      if(met.eq.0) irka(k,3)=1
      if(met.eq.1) then; krk1=k; itype1 = 'pqpq'; end if
      krk=k; itype = 'pqpq'
      call TimerStop('mrk_pqpq: generate array')

      call TimerStop('mrk_pqpq')

      END Subroutine mrk_pqpq


!======================================================================
      Subroutine diag_pqpq(k)
!======================================================================
!     Controls the scaling propeties for diagonal B-spline Rk-interals
!     (not implimented yet)
!----------------------------------------------------------------------
      Use DBS_grid, only: nv

      Implicit none
      Integer :: k,iv

      Do iv=1,nv; Call triang_pqpq(k,iv); End do

      END Subroutine diag_pqpq


!======================================================================
      Subroutine triang_pqpq(k,iv)
!======================================================================
!     Returns the two-dimensional array of B-spline integrals
!               <A_i B_j|r^k/r^(k+1)|C_i' D_j'>
!     over the given triangle diagonal cell
!
!     On entry   iv  -  index of the diagonal cell
!     --------
!
!     On exit    rkd(.,.,iv) - arrays of Rk B-spline integrals for given
!     --------                 interval iv in the reduced-dimension mode
!
!     Calls:   gauleg, zbsplvd
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_moments
      Use Timer

      Implicit none
      Integer :: iv,i,j, ip,jp, ii,jj, m, left, k
      Real(8) :: xbase
      Real(8) :: x(ks),w(ks), bp(ks),bq(ks)
      Real(8) :: bspTmp(ks,ksp)
      Real(8) :: bspTmq(ks,ksq)
      Real(8) :: INTP(ksp,ksp,ks)
      Real(8) :: INTQ(ksq,ksq,ks)

! ... setup the gaussian points

      Call gauleg(0.d0,1.d0,x,w,ks)

! ... first integration:

      left=iv+ks-1;  xbase=t(left)

      Call TimerStart("mrk_pqpq: First integration")
      DO m=1,ks

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks
       Call zbsplvd(nsp,ksp,nv, tp,iv+ksp-1, 1,gx(i),1,dbip)
       bspTmp(i,1:ksp)= dbip(1,1:ksp,1)
      END DO

      DO i=1,ks
       Call zbsplvd(nsq,ksq,nv, tq,iv+ksq-1, 1,gx(i),1,dbiq)
       bspTmq(i,1:ksq)= dbiq(1,1:ksq,1)
      END DO

! ... and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:)
      gx(:) = gw(:)*gx(:)**k

!            / r(iv,m)                             k
! ... INT =  |      bsp(iv,:,j)(r) bsp(iv,:,jp)(r) r  dr
!           / r_iv

      Do j=1,ksp
       gw(:) = gx(:)*bspTmp(:,j)
       Do jp=j,ksp; INTP(j,jp,m)= SUM(gw(:)*bspTmp(:,jp));  END DO
      End do

      Do j=1,ksq
       gw(:) = gx(:)*bspTmq(:,j)
       Do jp=j,ksq; INTQ(j,jp,m)= SUM(gw(:)*bspTmq(:,jp));  END DO
      End do

      END DO !  over m
      Call TimerStop("mrk_pqpq: First integration")

! ... second integration
      Call TimerStart("mrk_pqpq: Second integration")

      gx(:) = grw(iv,:)*grm(iv,:)**(k+1)

             rkd(:,:,iv) = 0.d0

      ii=0;  DO i=1,ksp;  DO ip=i,ksp; ii = ii+1
              bp(:) = pbsp(iv,:,i)*pbsp(iv,:,ip)*gx(:)
      jj=0;  DO j=1,ksq;  DO jp=j,ksq; jj = jj+1
              rkd(ii,jj,iv) = rkd(ii,jj,iv) + SUM(bp(:)*INTQ(j,jp,:))
             END DO; END DO
             END DO; END DO

      ii=0;  DO i=1,ksq;  DO ip=i,ksq; ii = ii+1
              bq(:) = qbsp(iv,:,i)*qbsp(iv,:,ip)*gx(:)
      jj=0;  DO j=1,ksp;  DO jp=j,ksp; jj = jj+1
              rkd(jj,ii,iv) = rkd(jj,ii,iv) + SUM(bq(:)*INTP(j,jp,:))
             END DO; END DO
             END DO; END DO
      Call TimerStop("mrk_pqpq: Second integration")

      End Subroutine triang_pqpq
