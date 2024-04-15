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
