!----------------------------------------------------------------------
!  Collection of simple sorting subroutines:
!
!  ISORT (n,S)       -  sorting an integer array S(1:n)
!  ISORT2(n,N1,N2)   -  sorting two array simultaneously
!  RSORT (n,R)       -  sorting real array R(1:n)
!  SORTA(n,S,IPT)    -  sorting the absolute value of real array S
!                       IPT - sorting pointer
!  SORTI(n,S,IPT)    -  sorting an integer array S(1:n)
!  SORT_IT(n,IS,IPT) -  sorting an integer array S(1:n) (more effective)
!  SORTR(n,S,IPT)    -  sorting real array S(1:n)
!  SORTI2(n,IS1,iS2,IPT) - sorting two array simultaneously
!----------------------------------------------------------------------

!======================================================================
      Integer Function ISORT (n,S)
!======================================================================
!     simple sorting for integer array NN(1:n)
!     isort - number of needed permutations
!----------------------------------------------------------------------
      Implicit none
      Integer :: n, i,j,k
      Integer :: S(n)
      ISORT = 0
      Do i=1,n-1
       Do j=i+1,n
        if(S(i).le.S(j)) Cycle
        k=S(i); S(i)=S(j); S(j)=k; ISORT=ISORT+1
       End do
      End do
      End Function ISORT


!======================================================================
      Integer Function ISORT2(n,N1,N2)
!======================================================================
!     sort two arrays simultaniously (fist N1, then N2)
!     and returns also the number of required permutations
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n
      Integer, intent(inout) :: N1(n),N2(n)
      Integer :: kz, i1,i2, k
      kz=0
      Do i1=1,n-1
       Do i2=i1+1,n
       if(n1(i1).gt.n1(i2).or.&
         (n1(i1).eq.n1(i2).and.n2(i1).gt.n2(i2)) ) then
        k=n1(i1); n1(i1)=n1(i2); n1(i2)=k
        k=n2(i1); n2(i1)=n2(i2); n2(i2)=k
        kz=kz+1
       end if
       End do
      End do
      ISORT2=kz
      End Function ISORT2


!======================================================================
      Integer Function RSORT (n,r)
!======================================================================
!     simple sorting real array r(1:n)
!     rsort - number of needed permutations
!----------------------------------------------------------------------
      Implicit none
      Integer :: n,i,j
      Real(8) :: r(n), rr
      RSORT = 0
      Do i=1,n-1
       Do j=i+1,n
        if(r(i).le.r(j)) Cycle
        rr=r(i); r(i)=r(j); r(j)=rr; RSORT=RSORT+1
       End do
      End do
      End Function RSORT


!====================================================================
      Subroutine SORTA(n,S,IPT)
!--------------------------------------------------------------------
!     provides sorting pointer IPT for absolute value of real array S
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in)  :: n
      Real(8), intent(in)  :: S(n)
      Integer, intent(out) :: IPT(n)
      Integer :: i,i1,j1,i2,j2
      Do i=1,n; IPT(i)=i;  End do
      Do i1=1,n-1; j1=IPT(i1)
       Do i2=i1+1,n;  j2=IPT(i2)
        if(abs(S(j1)).lt.abs(S(j2))) then
         IPT(i2)=j1; j1=j2
        end if
       End do
       IPT(i1)=j1
      End do
      End Subroutine SORTA

!====================================================================
      Subroutine SORTI(n,S,IPT)
!--------------------------------------------------------------------
!     gives IPT - sorting pointer for integer array S
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n
      Integer, intent(in), dimension(n) :: S
      Integer, intent(out), dimension(n) :: IPT
      Integer :: i,i1,j1,i2,j2

      Do i=1,n;  IPT(i)=i;  End do
      Do i1=1,n-1;  j1=IPT(i1)
       Do i2=i1+1,n; j2=IPT(i2)
        if(S(j1).gt.S(j2)) then
         IPT(i2)=j1; j1=j2
        end if
       End do
       IPT(i1)=j1
      End do
      End Subroutine SORTI


!====================================================================
      Subroutine SORT_IT(n,IS,IPT)
!--------------------------------------------------------------------
!     gives IPT - pointer on sorting of integer array IS
!     more effective for big arrays ???
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n, IS(n)
      Integer, intent(out) :: IPT(n)
      Integer :: i,k,l,m,it,mt

      IPT(1)=1
      Do i = 2,n; it=IS(i); k=1; l=i-1

       if(it.ge.IS(ipt(l))) then; k=i; go to 2; end if
       if(it.le.IS(ipt(k))) go to 2

    1 if(k.gt.l) go to 2
      m=(k+l)/2; mt=IS(IPT(m))
      if    (it.lt.mt) then;  l = m - 1
      elseif(it.gt.mt) then;  k = m + 1
      else;  k=m+1; go to 2
      end if
      go to 1
    2 Continue

       if(k.lt.i) then
        Do m = i,k+1,-1; IPT(m) = IPT(m-1); End do
       end if
       IPT(k) = i

      End do

      End Subroutine SORT_it


!====================================================================
      Subroutine SORTI2(n,IS1,iS2,IPT)
!--------------------------------------------------------------------
!     gives IPT - pointer on sorting of integer arrays IS1, IS2
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in ) :: n
      Integer, intent(in ) :: IS1(n),IS2(n)
      Integer, intent(out) :: IPT(n)
      Integer :: i,i1,j1,i2,j2

      Do i=1,n; IPT(i)=i; End do

      Do i1=1,n-1
       j1=IPT(i1)
       Do i2=i1+1,n
        j2=IPT(i2)
        if(IS1(j1).gt.IS1(j2)) then
         IPT(i2)=j1
         j1=j2
        elseif(IS1(j1).eq.IS1(j2)) then
         if(IS2(j1).gt.IS2(j2)) then
          IPT(i2)=j1
          j1=j2
         end if
        end if
       End do
       IPT(i1)=j1
      End do

      End Subroutine SORTI2

!====================================================================
      Subroutine SORTR(n,S,IPT)
!--------------------------------------------------------------------
!     provide sorting pointer IPT for real array S
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in ) :: n
      Real(8), intent(in ) :: S(n)
      Integer, intent(out) :: IPT(n)
      Integer :: i,i1,j1,i2,j2
      Do i=1,n; IPT(i)=i;  End do
      Do i1=1,n-1;    j1=IPT(i1)
       Do i2=i1+1,n;  j2=IPT(i2)
        if(S(j1).gt.S(j2)) then; IPT(i2)=j1; j1=j2; end if
       End do
       IPT(i1)=j1
      End do
      End Subroutine SORTR

