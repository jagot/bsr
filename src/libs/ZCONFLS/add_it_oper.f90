!=====================================================================
      Subroutine Add_it_oper(is,js)
!=====================================================================
! ... record in IT_oper what has been done for the given case: is,js 
!---------------------------------------------------------------------
      Use conf_LS,      only: JT_oper, noper
      Use symt_list_LS, only: IT_oper, ij 
      Use term_exp,     only: kt1,kt2, IP_kt1,IP_kt2

      Implicit none
      Integer :: is,js, it,jt,ik,jk, i, k
      Integer(8), external :: DEF_ij8

       k = 0
       Do ik=1,kt1;  it=IP_kt1(ik)
       Do jk=1,kt2;  jt=IP_kt2(jk)
        if(is.eq.js.and.it.gt.jt) Cycle;  k = k + 1
        ij=DEF_ij8(it,jt)
        Do i=1,noper; if(JT_oper(k,i).gt.0) IT_oper(i,ij)=1; End do
       End do; End do

      End Subroutine Add_it_oper
