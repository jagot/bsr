!---------------------------------------------------------------------
      Subroutine prepare_iort_jj(nu)
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
      Integer, intent(in) :: nu
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

! ... additional orthogonality conditions 

      if(nu.gt.0) Call read_orth_jj(nu)
      
      End Subroutine prepare_iort_jj


!======================================================================
      Subroutine read_orth_jj(nuc)
!======================================================================
!     read from file nuc and store in array IORT (module orb_jj)
!     the imposed orthogonal conditions given as 
!     < n1k1| n2k2>=0  (or 1 | 2)
!----------------------------------------------------------------------
      Use orb_jj

      Character(15) :: Aort

      rewind(nuc)
    1 read(nuc,'(a)',end=2) Aort
      if(Aort(1:1).ne.'<') go to 1 
      Call EL_nljk(Aort(2: 6),n1,kap1,l1,j1,k1)
      i1 = Ifind_jjorb(n1,kap1,k1,0)
      if(i1.eq.0.and.k1.ne.0) go to 1
      Call EL_nljk(Aort(8:12),n2,kap2,l2,j2,k2)
      i2 = Ifind_jjorb(n2,kap2,k2,0)
      if(i2.eq.0.and.k2.ne.0) go to 1
      if(kap1.ne.kap2) go to 1 
      read(Aort(15:15),'(i1)') ii
      if(k1.gt.0.and.k2.gt.0) then
       Call Iadd_orth(i1,i2,ii)
      elseif(k1.gt.0.and.k2.eq.0) then
       Do i2=1,nwf
        if(kap2.ne.KEF(i2).or.n2.ne.NEF(i2)) Cycle
        Call Iadd_orth(i1,i2,ii)
       End do       
      elseif(k1.eq.0.and.k2.gt.0) then
       Do i1=1,nwf
        if(kap1.ne.KEF(i1).or.n1.ne.NEF(i1)) Cycle
        Call Iadd_orth(i1,i2,ii)
      End do       
      elseif(k1.eq.0.and.k2.eq.0) then
       Do i1=1,nwf
        if(kap1.ne.KEF(i1).or.n1.ne.NEF(i1)) Cycle
       Do i2=1,nwf
        if(kap2.ne.KEF(i2).or.n2.ne.NEF(i2)) Cycle
        Call Iadd_orth(i1,i2,ii)
       End do; End do       
      end if

      go to 1
    2 Continue

      End Subroutine read_orth_jj 

