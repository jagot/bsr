!====================================================================
    Subroutine vinty_pq (ksm,ygr,yv,bsp)
!====================================================================
!   Computes the vector elements   <B_i, y(r)>
!--------------------------------------------------------------------
!   on entry
!   --------
!       ygr   array of values of a specific function  y(r) at the
!             gaussian points of each interval, weighted by the
!             gaussian weight
!
!   on exit
!   -------
!       yv    vector of integrals of <B_i, y(r)>, where i=1,..,ns
!
!--------------------------------------------------------------------
    Use DBS_grid

    Implicit none
    Integer, intent(in) :: ksm
    Real(8), intent(in)  :: ygr(nv,ks), bsp(nv+1,ks,ks)
    Real(8), intent(out) :: yv(nv+ksm-1)
    Integer :: iv, ith, i

    yv = 0.d0
    Do iv = 1,nv; Do ith = 1,ksm; i = iv+ith-1
      yv(i) = yv(i) + SUM(ygr(iv,:)*bsp(iv,:,ith))
    End Do; End Do

    End Subroutine vinty_pq


