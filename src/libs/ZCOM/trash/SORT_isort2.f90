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


