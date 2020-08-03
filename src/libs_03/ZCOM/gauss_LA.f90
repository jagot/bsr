!----------------------------------------------------------------------
      Subroutine gaussj(a,n,np,b,m,mp)
!----------------------------------------------------------------------
!     Linear equation solution by Gauss-Jordan elimination.
!     On input:
!     a(1:n,1:n) is an input matrix stored in an array of physical
!                dimensions np by np.
!     b(1:n,1:m) is an input matrix containing the m right-hand side
!                vectors, stored in an array np by mp.
!     On output:
!     a(1:n,1:n) is replaced by its matrix inverse
!     b(1:n,1:m) is replaced by the set of solution vectors.
!     (from numerical recipe,  www.nr.com)
!----------------------------------------------------------------------
      Integer :: m,mp,n,np
      Real(8) :: a(np,np),b(np,mp)
      Real(8) :: big,dum,pivinv, zero=0.d0, one=1.d0 
      Integer :: i,icol,irow,j,k,l,ll
      Integer :: indxc(n),indxr(n),ipiv(n)

! ... The integer arrays ipiv, indxr, and indxc are used
! ... for bookkeeping on the pivoting.

      ipiv = 0

! ... This is the main loop over the columns to be reduced.

      Do i=1,n;   big=zero

! ... This is the outer loop of the search for a pivot element.

      do j=1,n;  if(ipiv(j).eq.1) Cycle 
      do k=1,n;  if(ipiv(k).ne.0) Cycle
         if (abs(a(j,k)).ge.big) then
          big=abs(a(j,k)); irow=j; icol=k
         end if
      end do 
      end do 
      ipiv(icol)=ipiv(icol)+1

! ... We now have the pivot element, so we interchange rows, if needed, 
! ... to put the pivot element on the diagonal. The columns are not 
! ... physically interchanged, only relabeled:
! ... indxc(i), the column of the ith pivot element, is the ith column
! ...           that is reduced, while
! ... indxr(i)  is the row in which that pivot element was originally 
! ...           located. If indxr(i) /= indxc(i) there is an implied
! ... column interchange. With this form of bookkeeping, the solution
! ... b’s will end up in the correct order, and the inverse matrix will
! ... be scrambled by columns.

      if (irow.ne.icol) then
       do l=1,n
        dum=a(irow,l); a(irow,l)=a(icol,l); a(icol,l)=dum
       end do
       do l=1,m
        dum=b(irow,l); b(irow,l)=b(icol,l); b(icol,l)=dum
       end do 
      end if

! ... We are now ready to divide the pivot row by the pivot element,
! ... located at irow and icol.

      indxr(i)=irow 
      indxc(i)=icol 
      if (a(icol,icol).eq.zero) Stop 'singular matrix in gaussj'
      pivinv=one/a(icol,icol)
      a(icol,icol)=one
      do l=1,n; a(icol,l)=a(icol,l)*pivinv; end do
      do l=1,m; b(icol,l)=b(icol,l)*pivinv; end do

! ... Next, we reduce the rows, except for the pivot one, of course.

      do ll=1,n; if(ll.eq.icol) Cycle 
       dum=a(ll,icol)
       a(ll,icol)=zero 
       do l=1,n; a(ll,l)=a(ll,l)-a(icol,l)*dum; end do
       do l=1,m; b(ll,l)=b(ll,l)-b(icol,l)*dum; end do
      end do

     End do

! ... This is the end of the main loop over columns of the reduction.
! ... It only remains to unscramble the solution in view of the column
! ... interchanges. We do this by interchanging pairs of columns in 
! ... the reverse order that the permutation was built up.
   
     do l=n,1,-1; if(indxr(l).eq.indxc(l)) Cycle
     do k=1,n
      dum=a(k,indxr(l)); a(k,indxr(l))=a(k,indxc(l)); a(k,indxc(l))=dum
     end do 
     end do

     End Subroutine gaussj
