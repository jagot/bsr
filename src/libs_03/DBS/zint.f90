
!======================================================================
      Real(8) Function zint(int,i1,j1,i2,j2,k)
!======================================================================
!     selects different types of two-electron integrals      
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: int,k,i1,j1,i2,j2
      Real(8), external :: rk_pq, rk_pppp, rk_qqqq, rk_pqpq, rk_qpqp, &
                           sk_ppqq_pq, sk_pqqp_pq
      Select case(int)
       case(1);         zint = rk_pq(i1,j1,i2,j2,k)
       case(2);         zint = sk_ppqq_pq(i1,j1,i2,j2,k)
       case(3);         zint = rk_pppp(i1,j1,i2,j2,k)
       case(4);         zint = rk_qqqq(i1,j1,i2,j2,k)
       case(5);         zint = rk_pqpq(i1,j1,i2,j2,k)
       case(6);         zint = rk_qpqp(i1,j1,i2,j2,k)
       case(7);         zint = sk_ppqq_pq(i1,j1,i2,j2,k)
       case(8);         zint = sk_pqqp_pq(i1,j1,i2,j2,k)
       case default;    Stop 'zint: unknown int'
      End Select

      End Function zint

!======================================================================
      Real(8) Function zint_pq(int,i1,j1,i2,j2,k)
!======================================================================
!     selects different types of two-electron integrals      
!----------------------------------------------------------------------
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: int,k,i1,j1,i2,j2
      Real(8), external :: rk_pq, rk_pppp, rk_qqqq, rk_pqpq, rk_qpqp, &
                           sk_ppqq_pq, sk_pqqp_pq, Vp_dhl
      Select case(int)
       case(-1);        zint_pq = rk_pq(i1,j1,i2,j2,k)
       case(0);         zint_pq = Vp_dhl(kbs(i1),pq(1,1,i1),kbs(i2),pq(1,1,i2))
       case(1);         zint_pq = rk_pppp(i1,j1,i2,j2,k)
       case(2);         zint_pq = rk_qqqq(i1,j1,i2,j2,k)
       case(3);         zint_pq = rk_pqpq(i1,j1,i2,j2,k)
       case(4);         zint_pq = rk_qpqp(i1,j1,i2,j2,k)
       case(5);         zint_pq = sk_ppqq_pq(i1,j1,i2,j2,k)
       case(6);         zint_pq = sk_pqqp_pq(i1,j1,i2,j2,k)
       case default;    Stop 'zint_pq: unknown int'
      End Select

      End Function zint_pq
