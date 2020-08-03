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

      END SUBROUTINE SORT_it

