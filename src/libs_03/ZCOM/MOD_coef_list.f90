!====================================================================
      Module coef_list
!====================================================================
!     contains a set of coefficients ordering according to
!     "nnk" identifiers 
!     coefficients less than eps_C are ignored
!--------------------------------------------------------------------
      Implicit none 

      Integer :: nrk = 0       ! current number of coefficients
      Integer :: mrk = 0       ! maximum dimension
      Integer :: irk = 2**10   ! initial dimension
      Integer :: nnk = 0       ! number of identifiers          

! ... coefficients:

      Real(8), allocatable :: crk(:)   
    
! ... their attributes:

      Integer, allocatable :: krk(:,:)

      Real(8) :: eps_C = 1.d-12        
	  
      End MODULE coef_list


!======================================================================
      Subroutine alloc_coef_list(m,k,eps)
!======================================================================
      Use coef_list

      Implicit none
      Integer, intent(in) :: m,k
      Real(8), intent(in) :: eps
      Integer, allocatable :: ia(:,:)
      Real(8), allocatable :: ra(:)

      if(eps.gt.0.d0) eps_C=eps

      if(m.le.0.or.k.le.0) then
       if(allocated(crk)) Deallocate (crk,krk)
       nrk = 0; mrk = 0; nnk = 0 
      elseif(.not.allocated(crk)) then
       mrk = m; nrk = 0; nnk = k
       Allocate(crk(mrk),krk(k,m))
      elseif(m.le.mrk) then
       Return
      elseif(nrk.eq.0) then
       Deallocate (crk,krk)
       mrk = m; nnk = k
       Allocate(crk(mrk),krk(k,mrk))
      else
       if(k.gt.nnk) Stop 'alloc_coef_list: k > nnk'
       Allocate(ia(nnk,nrk))
       ia(:,1:nrk)=krk(:,1:nrk); Deallocate(krk)
       Allocate(krk(nnk,m)); krk(:,1:nrk)=ia(:,1:nrk)
       Deallocate(ia)
       Allocate(ra(nrk))
       ra(1:nrk)=crk(1:nrk); Deallocate(crk)
       Allocate(crk(m)); crk(1:nrk)=ra(1:nrk)
       Deallocate(ra)
       mrk = m
       write(*,'(a,2i10/a,f10.2,a)') ' Realloc_coef_list: new dimension = ', m,k, &
        ' memory requred = ', 4.d0*m*(2+nnk)/(1024*1024),'  Mb'
      end if

      End Subroutine alloc_coef_list


!======================================================================
      Subroutine Add_coef_list(n,in,C)
!======================================================================
!     add new data to the list 
!----------------------------------------------------------------------
      Use coef_list

      Implicit none
      Integer, intent(in) :: n,in(n)
      Real(8), intent(in) :: C
      Integer :: i,j,k,l,m

      if(abs(C).lt.eps_C) Return

      if(mrk.eq.0) Call alloc_coef_list(irk,n,0.d0)

! ... search position (k) for new integral

      k=1; l=nrk
    1 if(k.gt.l) go to 2              
      m=(k+l)/2
      
      j = 0
      Do i=1,n
       if(in(i).lt.krk(i,m)) then; l=m-1; Exit; end if
       if(in(i).gt.krk(i,m)) then; k=m+1; Exit; end if
       j = i
      End do

      if(j.eq.n) then    ! the same integral 
        crk(m)=crk(m)+C
        Return 
      end if

      go to 1
    2 Continue 

! ... shift the rest data up:

      Do i=nrk,k,-1
       m = i + 1;  crk(m) = crk(i);   krk(:,m) = krk(:,i)
      End do

! ... add new integral:

      crk(k)=C; krk(:,k)=in(:); nrk=nrk+1
      if(nrk.eq.mrk) Call alloc_coef_list(mrk+irk,nnk,0.d0) 

      End Subroutine Add_coef_list


