!----------------------------------------------------------------------
      Subroutine Sum_Term
!----------------------------------------------------------------------
!     exhaustion of shell-terms
!----------------------------------------------------------------------
      Use conf_LS, only: msh,no,ln,iq,LS

      Implicit none 
      Integer :: mt(msh),nt(msh)
      Integer :: i,i1,i2,ii,IA,IL,IS
      Integer, external :: Iterm_LS

!     mt(i) - the number of term in shell i
!     nt(i) - the term inder consideration

      i1=1                     ! i1 - low  limit of shells
      i2=no                    ! i2 - high limit of shells
      Do i=i1,i2
       mt(i)=Iterm_LS(ln(i),iq(i),-1,IA,IL,IS)
      End do

      i=i1                     ! first shell under consideration
      nt(i)=1

    1 ii=Iterm_LS(ln(i),iq(i),nt(i),LS(i,1),LS(i,2),LS(i,3))
      if(i.lt.i2) then
         i=i+1; nt(i)=1; go to 1
      else
         CALL Sum_Iterm
      end if

    2 nt(i)=nt(i)+1
      if(nt(i).gt.mt(i)) then
        if(i.eq.i1) go to 3
        i=i-1; go to 2
        end if
      go to 1

    3 Return

      End  Subroutine Sum_Term



!======================================================================
      Subroutine Sum_Iterm
!======================================================================
!     exhaustion of intermediate terms
!----------------------------------------------------------------------
      Use conf_LS, only: msh,no,LS

      Implicit none
      Integer :: LL_min(msh),LL_max(msh),SS_min(msh),SS_max(msh)
      Integer :: i,i1,i2,j1,j2

      LS(1,4)=LS(1,2)
      LS(1,5)=LS(1,3)
      if(no.eq.1) then;  CALL Save_term; Return; end if

      i1=2                         ! i1 - low  limit
      i2=no                        ! i2 - high limit in array LS(...)
      i=i1
    1 j1=i-1; j2=i

      LL_min(i)=IABS(LS(j1,4)-LS(j2,2))+1
      LL_max(i)=     LS(j1,4)+LS(j2,2) -1
      SS_min(i)=IABS(LS(j1,5)-LS(j2,3))+1
      SS_max(i)=     LS(j1,5)+LS(j2,3) -1
      LS(i,4)=LL_min(i)
      LS(i,5)=SS_min(i)

    2 if(i.lt.i2) then
       i=i+1; go to 1
      else
       CALL Save_term
      end if

    3 if(LS(i,5).lt.SS_max(i)) then
         LS(i,5)=LS(i,5)+2
         go to 2
      elseif(LS(i,4).lt.LL_max(i)) then
         LS(i,4)=LS(i,4)+2
         LS(i,5)=SS_min(i)
         go to 2
      else
         if(i.le.i1) go to 4
         i=i-1; go to 3
      end if

    4 Return

      End  Subroutine Sum_Iterm


!======================================================================
      Subroutine Save_term
!======================================================================
      Use conf_LS
       
      Implicit none
      Integer :: i, ILT,IST, j1,j2
      Integer, external :: Jadd_symt_LS

! ... check the total term: 

      ILT=LS(no,4)
      IST=LS(no,5)
      j1=iabs(ILT-IST)+1
      j2=iabs(ILT+IST)-1

      if(L_min.gt.0.and.L_min.gt.ILT) Return
      if(L_max.gt.0.and.L_max.lt.ILT) Return
      if(S_min.gt.0.and.S_min.gt.IST) Return
      if(S_max.gt.0.and.S_max.lt.IST) Return
      if(J_min.gt.0.and.J_min.gt. j2) Return
      if(J_max.gt.0.and.J_max.lt. j1) Return

! ... add configuration: 

      i = Jadd_symt_LS(iconf,no,LS)

      End  Subroutine Save_term


