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
