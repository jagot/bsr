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

