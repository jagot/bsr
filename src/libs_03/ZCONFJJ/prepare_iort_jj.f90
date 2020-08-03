!---------------------------------------------------------------------
      Subroutine prepare_iort_jj
!---------------------------------------------------------------------
!     prepares the orthogonality conditions (IORT) for all orbitals
!     according to the orbital set numbers KEF, and additional
!     conditions from c-file (unit nu)
!     iort(i,j) = 0    -  for orthogonal orbitals i,j
!     iort(i,j) = 2    -  for non-orthogonal orbitals i,j
!     iort(i,j) = 1    -  for the normalized orbitals 
!     OPTIONS:
!     JORT<0 - full orthogonality
!     JORT=0 - full non-orthogonality
!     JORT>0 - partial orthogonality, i.e.
!     the orbitals  are orthogonal within one set with the same KEF,
!     and all orbitals are orthonomal to the orbitals with KEF=0
!---------------------------------------------------------------------
      Use orb_jj

      Implicit none
      Integer :: i1,i2
      
      Do i1=1,nwf
      Do i2=1,i1
       if(KEF(i1).ne.KEF(i2)) Cycle
       if(JORT.lt.0) then             ! full orthogonality
        if(i1.eq.i2) Call Iadd_orth(i1,i2,1)
       elseif(JORT.eq.0) then         ! full non-orthogonality
        Call Iadd_orth(i1,i2,2) 
       elseif(JORT.gt.0) then         ! partial orthogonality
        if(IEF(i1).eq.IEF(i2)) then
         if(NEF(i1).eq.NEF(i2)) Call Iadd_orth(i1,i2,2)          
        elseif(IEF(i1)*IEF(i2).eq.0) then
         Cycle
        else
         Call Iadd_orth(i1,i2,2) 
        end if
       end if
      End do
      End do

      End Subroutine prepare_iort_jj


