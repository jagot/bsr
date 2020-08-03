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
