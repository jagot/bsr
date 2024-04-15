!======================================================================
      Module tm_LS
!======================================================================
!     Containes the allocatable set of records with three real and five
!     integer  parameters
!     Used to save T-matrix elements (for one transition):
!     STpar  -  total spin of the given partial wave
!     LTpar  -  total orbital momemtum of the given partial wave
!     l1, l2 -  l-values for ingoing and outgoing electron
!     ip     -
!     EK     -  electron energy
!     TR, TI -  real and imaginary parts of T-matrix element
!     jp     -
!     kcase  -
!----------------------------------------------------------------------
      Implicit none

      Integer :: ndata = 0       ! number of entries
      Integer :: mdata = 0       ! current dimentsion of list
      Integer :: idata = 2**10   ! supposed max. dimentsion

      Integer, allocatable :: STpar(:),LTpar(:),l1(:),l2(:),ip(:)

      Real(8), allocatable :: EK(:),TR(:),TI(:)

      Integer, allocatable :: jp(:,:)

      Integer :: kcase

      End Module tm_LS


!======================================================================
      Subroutine Alloc_tm_LS(m)
!======================================================================
!     allocate or deallocate the data in "tm_LS"
!----------------------------------------------------------------------
      Use tm_LS

      Implicit none
      Integer :: m,n
      Integer, allocatable :: iarr(:)
      Real(8), allocatable :: rarr(:)

      if(m.le.0) then
       if(allocated(EK)) Deallocate(EK,TR,TI,LTpar,ip,l1,l2,STpar)
       mdata = 0; ndata = 0
      elseif(.not.allocated(EK)) then
       Allocate(EK(m),TR(m),TI(m),LTpar(m),ip(m),l1(m),l2(m),STpar(m))
       mdata = m; ndata = 0; ip=0
      elseif(m.le.mdata) then
       Return
      elseif(ndata.gt.0) then
       n = ndata
       Allocate(iarr(n))
       iarr(1:n)=LTpar(1:n); Deallocate(LTpar); Allocate(LTpar(m))
       LTpar(1:n)=iarr(1:n)
       iarr(1:n)=ip(1:n); Deallocate(ip); Allocate(ip(m))
       ip(1:n)=iarr(1:n)
       iarr(1:n)=l1(1:n); Deallocate(l1); Allocate(l1(m))
       l1(1:n)=iarr(1:n)
       iarr(1:n)=l2(1:n); Deallocate(l2); Allocate(l2(m))
       l2(1:n)=iarr(1:n)
       iarr(1:n)=STpar(1:n); Deallocate(STpar); Allocate(STpar(m))
       STpar(1:n)=iarr(1:n)
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

      End Subroutine alloc_tm_LS


!=======================================================================
      Subroutine Add_tm_LS(e,ar,ai,L,S,m1,m2,ii)
!=======================================================================
!     add (or substitude) data to the list
!-----------------------------------------------------------------------
      Use tm_LS

      Implicit none
      Integer :: L,S,m1,m2,ii, i
      Real(8) :: e,ar,ai

      if(mdata.eq.0) Call Alloc_tm_LS(idata)

! ... check if the same data are already in the list:

      Do i=1,ndata
       if(L.ne.LTpar(i)) Cycle
       if(S.ne.STpar(i)) Cycle
       if(ii.ne.ip(i)) Cycle
       if(m1.ne.l1(i)) Cycle
       if(m2.ne.l2(i)) Cycle
       if(e.ne.EK(i)) Cycle
       TR(i)=ar
       TI(i)=ai
       Return
      End do

! ... add new record:

      if(ndata.eq.mdata) Call Alloc_tm_LS(mdata+idata)

      ndata = ndata+1
      LTpar(ndata) = L
      STpar(ndata) = S
      ip(ndata) = ii
      l1(ndata) = m1
      l2(ndata) = m2
      EK(ndata) = e
      TR(ndata) = ar
      TI(ndata) = ai

      End Subroutine Add_tm_LS


!======================================================================
      Subroutine Sort_tm_LS
!======================================================================
!     sortin data in tm_LS
!----------------------------------------------------------------------
      Use tm_LS

      Implicit none
      Integer :: i,j

      Do i=1,ndata-1;  Do j=i+1,ndata

       if(EK(i).lt.EK(j)) Cycle
       if(EK(i).gt.EK(j)) then; Call  Change_tm_LS(i,j); Cycle; end if

       if(ip(i).lt.ip(j)) Cycle
       if(ip(i).gt.ip(j)) then; Call  Change_tm_LS(i,j); Cycle; end if

       if(STpar(i).lt.STpar(j)) Cycle
       if(STpar(i).gt.STpar(j)) then; Call  Change_tm_LS(i,j); Cycle; end if

       if(LTpar(i).lt.LTpar(j)) Cycle
       if(LTpar(i).gt.LTpar(j)) then; Call  Change_tm_LS(i,j); Cycle; end if

       if(l1(i).lt.l1(j)) Cycle
       if(l1(i).gt.l1(j)) then; Call  Change_tm_LS(i,j); Cycle; end if

       if(l2(i).lt.l2(j)) Cycle
       if(l2(i).gt.l2(j)) then; Call  Change_tm_LS(i,j); Cycle; end if

      End do; End do

      End Subroutine Sort_tm_LS


!=======================================================================
      Subroutine Change_tm_LS(i,j)
!=======================================================================
!     exchange two records "i" and "j" in tmat_lis
!------------------------------------------------------------------------
      Use tm_LS
      Implicit none
      Integer :: i,j,k
      Real(8) :: S

      S = EK(i); EK(i)=EK(j); EK(j)=S
      S = TR(i); TR(i)=TR(j); TR(j)=S
      S = TI(i); TI(i)=TI(j); TI(j)=S

      k = LTpar(i); LTpar(i)=LTpar(j); LTpar(j)=k
      k = STpar(i); STpar(i)=STpar(j); STpar(j)=k
      k = ip(i);   ip(i)=ip(j);     ip(j)=k
      k = l1(i);   l1(i)=l1(j);     l1(j)=k
      k = l2(i);   l2(i)=l2(j);     l2(j)=k

      End Subroutine  Change_tm_LS
