!=====================================================================
      Module term_LS
!=====================================================================
!     description of involved terms 
!     ROUTINES:  alloc_term_LS
!                ifind_term()
!---------------------------------------------------------------------
      Implicit none

      Integer :: nterms =  0      ! current number of terms
      Integer :: mterms =  0      ! max. number of terms
      Integer :: iterms = 10      ! initial suggestion for mterm

      Integer, allocatable :: ILterm(:)   ! total L
      Integer, allocatable :: ISterm(:)   ! total S

      Integer :: ncfgt = 0
      Integer, allocatable :: LSP(:)      ! term pointer 

! ... first copy:

      Integer :: nterm1 =  0      
      Integer, allocatable :: ILterm1(:)  
      Integer, allocatable :: ISterm1(:)  
      Integer :: ncfgt1 = 0
      Integer, allocatable :: LSP1(:)   

! ... second copy:

      Integer :: nterm2 =  0      
      Integer, allocatable :: ILterm2(:)  
      Integer, allocatable :: ISterm2(:)  
      Integer :: ncfgt2 = 0
      Integer, allocatable :: LSP2(:)   

! ... DJ-factors

      Real(8), allocatable :: DJ(:,:), DJM(:,:)

      End Module term_LS


!======================================================================
      Subroutine alloc_term_LS(m)
!======================================================================
!     allocate arrays in module term_LS
!----------------------------------------------------------------------
      Use term_LS

      Implicit none
      Integer :: m
      Integer, allocatable :: ia(:)

      if(m.le.0) then
       if(allocated(ILterm)) Deallocate (ILterm,ISterm)
       mterms = 0; nterms = 0
      elseif(.not.allocated(ILterm)) then
       mterms = m
       Allocate(ILterm(mterms),ISterm(mterms))
       nterms = 0
      elseif(m.le.mterms) then
       Return
      elseif(nterms.le.0) then
       Deallocate (ILterm,ISterm); nterms = 0
       mterms = m
       Allocate(ILterm(mterms),ISterm(mterms))
      else
       mterms=m
       Allocate(ia(nterms))
       ia = ILterm(1:nterms); Deallocate(ILterm)
       Allocate(ILterm(mterms)); ILterm(1:nterms)=ia
       ia = ISterm(1:nterms); Deallocate(ISterm)
       Allocate(ISterm(mterms)); ISterm(1:nterms)=ia
       Deallocate(ia)
      end if

      End Subroutine alloc_term_LS


!=======================================================================
      Integer FUNCTION Ifind_term(L,S,job)
!=======================================================================
!     find term LS in the list, = 0, if no term
!     job = 0  -  no further actions
!     job = 1  -  stop if fail to find
!     job = 2  -  add new term
!------------------------------------------------------------------------
      Use term_LS

      Implicit none
      Integer, intent(in) :: L,S,job
      Integer :: i

      Ifind_term=0

      Do i=1,nterms
       if(L.ne.ILterm(i)) Cycle
       if(S.ne.ISterm(i)) Cycle
       Ifind_term = i
       Return
      End do
      if(job.eq.0) Return

      if(job.eq.1) then
       write(*,'(a,a,5i5)') 'Ifind_term can not find the term:',&
                            ' L,S = ',L,S
       Stop 
      end if

      if(nterms.ge.mterms) Call Alloc_term_LS(mterms+iterms)
      nterms = nterms + 1
      ILterm(nterms) = L
      ISterm(nterms) = S
      Ifind_term = nterms

      End Function Ifind_term


!=======================================================================
      Subroutine Shift_term1
!=======================================================================
!     save term information in the first set
!-----------------------------------------------------------------------
      Use term_LS

      if(nterms.eq.0) Return
      if(allocated(ILterm1)) Deallocate(ILterm1,ISterm1)
      nterm1=nterms
      Allocate(ILterm1(nterm1),ISterm1(nterm1))
      ILterm1 = ILterm(1:nterm1)
      ISterm1 = ISterm(1:nterm1)

      if(ncfgt.eq.0) Return 
      if(allocated(LSP1)) Deallocate(LSP1)
      ncfgt1=ncfgt 
      Allocate(LSP1(ncfgt1))
      LSP1 = LSP(1:ncfgt1)

      End Subroutine Shift_term1


!=======================================================================
      Subroutine Shift_term2
!=======================================================================
!     save term information in the second set
!-----------------------------------------------------------------------
      Use term_LS

      if(nterms.eq.0) Return
      if(allocated(ILterm2)) Deallocate(ILterm2,ISterm2)
      nterm2=nterms
      Allocate(ILterm2(nterm2),ISterm2(nterm2))
      ILterm2 = ILterm(1:nterm2)
      ISterm2 = ISterm(1:nterm2)

      if(ncfgt.eq.0) Return 
      if(allocated(LSP2)) Deallocate(LSP2)
      ncfgt2=ncfgt 
      Allocate(LSP2(ncfgt2))
      LSP2 = LSP(1:ncfgt2)

      End Subroutine Shift_term2


!=======================================================================
      Subroutine Alloc_LSP(m)
!=======================================================================
!     allocate term pointers for the given set of configurations
!-----------------------------------------------------------------------
      Use term_LS

      Implicit none
      Integer :: m

      if(allocated(LSP)) Deallocate(LSP)
      ncfgt=m
      if(m.le.0) Return
      Allocate(LSP(ncfgt))

      End Subroutine Alloc_LSP


