!=======================================================================
      Subroutine INT_de (p1,p2,x,int,icase,sym_x)
!=======================================================================
!     Contribution to channel-channel interaction matrix from
!     the two-electron integrals:
!
!     icase = 1 -->  INT( . p1; . p2)
!     icase = 2 -->  INT( p1 .;p2 . )
!     icase = 3 -->  INT( . p1;p2 . )
!     icase = 4 -->  INT( p1 .; . p2)
!
!     result is in x array: x/ich,jch/
!-----------------------------------------------------------------------
      Use spline_param

      Implicit none

      Integer, intent(in) :: int, icase
      Real(8), intent(in) :: p1(ns),p2(ns)
      Real(8), intent(out) :: x(ns,*)
      Character(1), intent(out) :: sym_x

      Real(8) :: d(ns,ns)
      Character(1) :: sym_i, sym_j, sym_d

      if(int.eq.1) then              !  Uk
       sym_i = 's'
       sym_j = 'n'
      elseif(int.eq.2) then          !  QK
       sym_i = 'n'
       sym_j = 'n'
      elseif(int.eq.3) then          !  Tk
       sym_i = 'n'
       sym_j = 'n'
      elseif(int.eq.4) then          !  Mk
       sym_i = 's'
       sym_j = 's'
      elseif(int.eq.5) then          !  Rk
       sym_i = 's'
       sym_j = 's'
      elseif(int.eq.8) then          !  Nk
       sym_i = 's'
       sym_j = 's'
      elseif(int.eq.9) then          !  Vk
       sym_i = 'n'
       sym_j = 's'
      elseif(int.eq.10) then         !  Nk
       sym_i = 's'
       sym_j = 's'
      else
       Stop ' INT_d:  unknown type of integral '
      end if

      sym_d = 'x'
      if(icase.eq.1) sym_d = sym_j
      if(icase.eq.2) sym_d = sym_i

      Call Density (ns,ks,d,p1,p2,sym_d)

      Call Convol  (ns,ks,x,d,icase,sym_i,sym_j)

      sym_x = 'x'
      if(icase.eq.1) sym_x = sym_i
      if(icase.eq.2) sym_x = sym_j

      End Subroutine INT_de

