!======================================================================
      Subroutine INT_v(i1,j1,i2,j2,k,int,icase,v)
!======================================================================
!                   k                    k
!     Evaluates  INT ( . j1; i2 j2),  INT ( i1 . ; i2 j2),
!                   k                    k
!                INT (i1 j1;  . j2),  INT ( i1 i2; i2  .).
!----------------------------------------------------------------------
      Use spline_param
      Use spline_orbitals, p => pbs

      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k,int,icase
      REAL(8), intent(out) :: v(ns)
      REAL(8)  :: x(ns,2*ks-1)
      Character(1) :: sym_x

      if(icase.eq.1) then
        CALL INT_de (p(1,j1),p(1,j2),x,int,1,sym_x)
        CALL bav (ns,ks,x,p(1,i2),v,sym_x,'r')
      elseif(icase.eq.2) then
        CALL INT_de (p(1,i1),p(1,i2),x,int,2,sym_x)
        CALL bav (ns,ks,x,p(1,j2),v,sym_x,'r')
      elseif(icase.eq.3) then
        CALL INT_de (p(1,j1),p(1,j2),x,int,1,sym_x)
        CALL bav (ns,ks,x,p(1,i1),v,sym_x,'l')
      elseif(icase.eq.4) then
        CALL INT_de (p(1,i1),p(1,i2),x,int,2,sym_x)
        CALL bav (ns,ks,x,p(1,j1),v,sym_x,'l')
      end if

      End Subroutine INT_v
