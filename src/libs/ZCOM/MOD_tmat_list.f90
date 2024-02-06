!======================================================================
      Module tmat_list
!======================================================================
!     Containes the allocatable data with three real and six integer
!     parameters
!----------------------------------------------------------------------
      Implicit none

      Integer :: ndata = 0       ! number of determinants
      Integer :: mdata = 0       ! current dimentsion of list
      Integer :: idata = 2**10   ! supposed max. dimentsion

      Integer, allocatable :: jjpar(:),ip(:),l1(:),l2(:),jk1(:),jk2(:)
      Real(8), allocatable :: EK(:),TR(:),TI(:)

      End Module tmat_list


!======================================================================
      Subroutine alloc_tmat_list(m)
!======================================================================
!     allocate or deallocate the data
!----------------------------------------------------------------------
      Use tmat_list

      Implicit none
      Integer :: m,n
      Integer, allocatable :: iarr(:)
      Real(8), allocatable :: rarr(:)

      if(m.le.0) then
       if(allocated(EK)) Deallocate(EK,TR,TI,jjpar,ip,l1,l2,jk1,jk2)
       mdata = 0; ndata = 0
      elseif(.not.allocated(EK)) then
       Allocate(EK(m),TR(m),TI(m),jjpar(m),ip(m),l1(m),l2(m),jk1(m),jk2(m))
       mdata = m; ndata = 0
      elseif(m.le.mdata) then
       Return
      elseif(ndata.gt.0) then
       n = ndata
       Allocate(iarr(n))
       iarr(1:n)=jjpar(1:n); Deallocate(jjpar); Allocate(jjpar(m))
       jjpar(1:n)=iarr(1:n)
       iarr(1:n)=ip(1:n); Deallocate(ip); Allocate(ip(m))
       ip(1:n)=iarr(1:n)
       iarr(1:n)=l1(1:n); Deallocate(l1); Allocate(l1(m))
       l1(1:n)=iarr(1:n)
       iarr(1:n)=l2(1:n); Deallocate(l2); Allocate(l2(m))
       l2(1:n)=iarr(1:n)
       iarr(1:n)=jk1(1:n); Deallocate(jk1); Allocate(jk1(m))
       jk1(1:n)=iarr(1:n)
       iarr(1:n)=jk2(1:n); Deallocate(jk2); Allocate(jk2(m))
       jk2(1:n)=iarr(1:n)
       Deallocate(iarr)
       Allocate(rarr(n))
       rarr(1:n)=EK(1:n); Deallocate(EK); Allocate(EK(m))
       EK(1:n)=rarr(1:n)
       rarr(1:n)=TR(1:n); Deallocate(TR); Allocate(TR(m))
       TR(1:n)=rarr(1:n)
       rarr(1:n)=TI(1:n); Deallocate(TI); Allocate(TI(m))
       TI(1:n)=rarr(1:n)
       Deallocate(rarr)
       mdata = m
      end if

      End Subroutine alloc_tmat_list


!=======================================================================
      Subroutine Add_tmat_list(e,ar,ai,jj,ii,m1,m2,k1,k2)
!=======================================================================
!     add (sum to) data to the list
!     we may need summation due to transition to other coupling
!-----------------------------------------------------------------------
      Use tmat_list

      Implicit none
      Integer :: jj,ii,m1,m2,k1,k2, i
      Real(8) :: e,ar,ai

      if(mdata.eq.0) Call Alloc_tmat_list(idata)

! ... check if the same data are already in the list:

      Do i=1,ndata
       if(jj.ne.jjpar(i)) Cycle
       if(ii.ne.ip(i)) Cycle
       if(m1.ne.l1(i)) Cycle
       if(m2.ne.l2(i)) Cycle
       if(k1.ne.jk1(i)) Cycle
       if(k2.ne.jk2(i)) Cycle
       if(e.ne.EK(i)) Cycle
       TR(i)=TR(i)+ar
       TI(i)=TI(i)+ai
       Return
      End do

! ... add new integral:

      if(ndata.eq.mdata) Call Alloc_tmat_list(mdata+idata)

      ndata = ndata+1
      jjpar(ndata) = jj
      ip(ndata) = ii
      l1(ndata) = m1
      l2(ndata) = m2
      jk1(ndata) = k1
      jk2(ndata) = k2
      EK(ndata) = e
      TR(ndata) = ar
      TI(ndata) = ai

      End Subroutine Add_tmat_list


!======================================================================
      Subroutine Sort_tmat_list
!======================================================================
!     sortin data in tmat_list
!----------------------------------------------------------------------
      Use tmat_list

      Implicit none
      Integer :: i,j

      Do i=1,ndata-1;  Do j=i+1,ndata

       if(EK(i).lt.EK(j)) Cycle
       if(EK(i).gt.EK(j)) then; Call Change_ij(i,j); Cycle; end if

       if(jjpar(i).lt.jjpar(j)) Cycle
       if(jjpar(i).gt.jjpar(j)) then; Call Change_ij(i,j); Cycle; end if

       if(ip(i).lt.ip(j)) Cycle
       if(ip(i).gt.ip(j)) then; Call Change_ij(i,j); Cycle; end if

       if(l1(i).lt.l1(j)) Cycle
       if(l1(i).gt.l1(j)) then; Call Change_ij(i,j); Cycle; end if

       if(jk1(i).lt.jk1(j)) Cycle
       if(jk1(i).gt.jk1(j)) then; Call Change_ij(i,j); Cycle; end if

       if(l2(i).lt.l2(j)) Cycle
       if(l2(i).gt.l2(j)) then; Call Change_ij(i,j); Cycle; end if

       if(jk2(i).lt.jk2(j)) Cycle
       if(jk2(i).gt.jk2(j)) then; Call Change_ij(i,j); Cycle; end if

      End do; End do


      End Subroutine Sort_tmat_list


!=======================================================================
      Subroutine Change_ij(i,j)
!=======================================================================
!     exchange two records "i" and "j" in tmat_lis
!------------------------------------------------------------------------
      Use tmat_list

      Implicit none
      Integer :: i,j,k
      Real(8) :: S

      S = EK(i); EK(i)=EK(j); EK(j)=S
      S = TR(i); TR(i)=TR(j); TR(j)=S
      S = TI(i); TI(i)=TI(j); TI(j)=S

      k = jjpar(i); jjpar(i)=jjpar(j); jjpar(j)=k
      k = ip(i);    ip(i)=ip(j);       ip(j)=k
      k = l1(i);    l1(i)=l1(j);       l1(j)=k
      k = l2(i);    l2(i)=l2(j);       l2(j)=k
      k = jk1(i);   jk1(i)=jk1(j);     jk1(j)=k
      k = jk2(i);   jk2(i)=jk2(j);     jk2(j)=k

      End Subroutine Change_ij
