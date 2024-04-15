!---------------------------------------------------------------------
      Subroutine Pre_iort(nu,nmax)
!---------------------------------------------------------------------
!     prepares the orthogonality conditions (IORT) for all orbitals
!     according to the orbital set numbers KEF, and additional
!     conditions from c-file (unit nu)
!
!     It is assumed that the orbitals of one set are orthogonal,
!     and all orbitals are orthonomal to the ones with KEF=0
!
!     JORT<0 - full orthogonality
!     JORT=0 - full non-orthogonality
!     JORT>0 - partial orthogonality
!
!     IORT=0 - orthogonal
!     IORT=1 - the same (overlap =1)
!     iort=2 - non-orthogonal
!---------------------------------------------------------------------
      Use orb_LS

      Implicit none
      Integer, intent(in) :: nu,nmax
      Integer :: i1,i2

      Do i1=1,nwf
      Do i2=1,i1

       IORT(i1,i2)=0
       if(LEF(i1).ne.LEF(i2)) Cycle

       if(JORT.lt.0) then             ! full orthogonality

        IORT(i1,i2)=0
        if(i1.eq.i2) IORT(i1,i2)=1

       elseif(JORT.eq.0) then         ! full non-orthogonality

        IORT(i1,i2)=2

       elseif(JORT.gt.0) then         ! partial orthogonality

        IORT(i1,i2)=2

        if(KEF(i1).ne.KEF(i2)) then
         if(KEF(i1)*KEF(i2).eq.0) IORT(i1,i2)=0
        else
         IORT(i1,i2)=0; if(NEF(i1).eq.NEF(i2)) IORT(i1,i2)=1
         if(nmax.gt.0) then
          if(nmax.gt.0.and.NEF(i1).gt.nmax) IORT(i1,i2)=2
          if(nmax.gt.0.and.NEF(i2).gt.nmax) IORT(i1,i2)=2
         end if
        end if

       end if

      End do
      End do

! ... additional orthogonality conditions

      if(nu.gt.0) Call R_orth(nu)   ! ???

      End Subroutine Pre_iort
