!====================================================================
      Integer Function ML_id (id)
!--------------------------------------------------------------------
!     gives orbital magnetic numbers for id'th orbital
!     in the (2j+1) representation, the following order being assumed:
!
!     0 -, 0 +, -1 -, -1 +, 1 -, 1 +, -2 -, -2 +, 2 -, 2 +, ...
!---------------------------------------------------------------------

      Integer :: k,ml

      if(id.lt.1) Stop ' ML_id: id<1'
      ml=(id-1)/4
      ml=ml+ml+1
      k=mod(id-1,4)
      if(k.eq.2.or.k.eq.3) ml=-ml
      ml_id=ml

      End Function ML_id


!--------------------------------------------------------------------
      Integer Function MS_id (id)
!--------------------------------------------------------------------
!     gives spin magnetic numbers for id'th orbital
!     in the (2j+1) representation, the following order being assumed:
!
!     0 -, 0 +, -1 -, -1 +, 1 -, 1 +, -2 -, -2 +, 2 -, 2 +, ...
!---------------------------------------------------------------------

      Integer, Intent(in) :: id

      if(id.lt.1) Stop ' MS_id: id<1'
      ms_id=mod(id-1,2)*2

      End Function MS_id
