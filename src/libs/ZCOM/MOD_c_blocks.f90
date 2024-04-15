!======================================================================
      Module c_blocks
!======================================================================
! ... contains a set of coefficients with four identifiers
! ... (k1,k2,k3,k4) and one pointer (ipt), devided in blocks;
! ... each block contain data for given type
!
! ... internal procedures: Alloc_c_block   (ntype,...)
!                          Add_coef_cblock (itype,C,j1,j2,j3,j4)
!                          Add_c_block     (ip,jp,C,j1,j2,j3,j4)
!                          Collect_coef    (itype)
!                          Merge_c_block   (nn, ip, jp, nc)
! ... external calls:      Do_coef         (itype)
!----------------------------------------------------------------------
      Implicit none

      Integer :: ncdata = 0             ! number of coefficients
      Real(8), allocatable :: CDATA(:)  ! coefficient values

! ... their attributes:

      Integer, allocatable :: K1(:),K2(:),K3(:),K4(:),IPT(:)

! ... default dimension limits:

      Integer :: nblock  =   1000         ! number of blocks in c_data
      Integer :: mblock  =   2000         ! size of blocks
      Integer :: kblock  =    100         ! max.nb for given type
      Integer :: ntype   =    100         ! number of integral types

      Real(8) :: mem_cdata  = 0.d0

! ... tolerence parameters:

      Real(8) :: eps_cdata  = 1.d-10      ! tolerance for coefficients

! ... pointer on first(last) element in given block

      Integer, allocatable :: ipblk(:), ipi(:)
      Integer, allocatable :: jpblk(:), ipj(:)

! ... current block for given type:

      Integer, allocatable :: iblk(:)

! ... number of blocks for given type:

      Integer, allocatable :: nblk(:)

! ... chain ponter for given type:

      Integer, allocatable :: kblk(:,:)

      End Module c_blocks


!======================================================================
      Subroutine alloc_c_blocks(nt,mb,nb,kb,eps_c)
!======================================================================
!     allocate (deallocate) arrays in module c_blocks
!----------------------------------------------------------------------
      Use c_blocks

      Implicit none
      Integer, intent(in) :: nt,mb,nb,kb
      Real(8), intent(in) :: eps_c
      Integer :: m,i,j,k

      if(allocated(CDATA)) Deallocate(CDATA,K1,K2,K3,K4,IPT, &
                           ipblk,jpblk,ipi,ipj,iblk,nblk,kblk)
      mem_cdata = 0.0
      if(nt.lt.0)  Return

      if(nt .gt. 0) ntype   =  nt
      if(mb .gt. 0) mblock  =  mb
      if(nb .gt. 0) nblock  =  nb
      if(kb .gt. 0) kblock  =  kb

      if(eps_c.gt.0.d0) eps_cdata = eps_c

      m = mblock*nblock
      Allocate(CDATA(m),K1(m),K2(m),K3(m),K4(m),IPT(m))
      Allocate(ipblk(nblock),jpblk(nblock),ipi(nblock),ipj(nblock), &
               iblk(ntype),nblk(ntype), kblk(ntype,nblock) )

      m = 7*m + 4*nblock + (2+nblock)*ntype
      mem_cdata = m * 4.0 / (1024 * 1024)

! ... initilize all blocks:

      Do i=1,nblock;  ipblk(i)=(i-1)*mblock+1; jpblk(i)=-1; End do

      CDATA = 0.d0

! ... assign one block to each type:

      i = 0
      Do j = 1,ntype
       i = i + 1
       iblk(j)=i; nblk(j)=1; kblk(j,1)=i; jpblk(i)=0
      End do

      End Subroutine alloc_c_blocks


!======================================================================
      Subroutine realloc_c_blocks(nt)
!======================================================================
!     re-allocate  arrays in module c_blocks
!----------------------------------------------------------------------
      Use c_blocks

      Implicit none
      Integer, intent(in) :: nt
      Integer, allocatable :: iar(:), jar(:,:)
      Integer :: i,m,itype,jtype

      if(nt.le.ntype) Return

      Allocate(iar(ntype))
      iar = iblk; Deallocate(iblk); Allocate(iblk(nt)); iblk=0
      iblk(1:ntype) = iar
      iar = nblk; Deallocate(nblk); Allocate(nblk(nt)); nblk=0
      nblk(1:ntype) = iar
      Deallocate(iar)

      Allocate(jar(ntype,nblock))
      jar = kblk; Deallocate(kblk); Allocate(kblk(nt,nblock)); kblk=0
      kblk(1:ntype,:) = jar
      Deallocate(jar)

! ... assign one block to new type:

      Do itype = ntype+1,nt

       m = 0
       Do i=1,nblock
        if(jpblk(i).ge.0) Cycle
        iblk(itype) = i
        nblk(itype) = 1
        kblk(itype,1) = i
        jpblk(i) = 0
        m = 1; Exit
       End do

       if(m.eq.1) Cycle

! ... everything is full - it is time to generate matrix;
! ... we do it for type with biggest occupation

       jtype = maxloc(nblk,DIM=1)
       Call Collect_coef(jtype)

! ... find again new block:

       m = 0
       Do i=1,nblock
        if(jpblk(i).ge.0) Cycle
        iblk(itype) = i
        nblk(itype) = 1
        kblk(itype,1) = i
        jpblk(i) = 0
        m = 1; Exit
       End do

       if(m.eq.0) Stop ' Realloc_c_bloks: problems with new block '

      End do

      ntype = nt
      m = 7*mblock*nblock + 4*nblock + (2+nblock)*ntype
      mem_cdata = m * 4.0 / (1024 * 1024)

      End Subroutine realloc_c_blocks


!======================================================================
      Subroutine Add_coef_cblock (C,j1,j2,j3,j4,itype)
!======================================================================
!     add new coefficient to the list for given type
!----------------------------------------------------------------------
      Use c_blocks

      Implicit none
      Real(8), intent(in) :: C
      Integer, intent(in) :: j1,j2,j3,j4,itype
      Integer :: i,k,m,n,ip,jp,jtype, m1,m2,m3

      if(itype.gt.ntype)  Call realloc_c_blocks(ntype+10)

! Call Decode_type(itype,m1,m2,m3)
! write(81,*) 'Add_coef_cblock',itype, m1,m2,m3

! ... add coefficient to list:

      i = iblk(itype); ip = ipblk(i); jp = jpblk(i)
      Call Add_c_block(ip,jp,C,j1,j2,j3,j4)
      jpblk(i) = jp

! ... check if the block full:

      if(jp-ip+1.lt.mblock) Return

! ... find new block:

      m = 0
      if(nblk(itype).lt.kblock) then
       Do i=1,nblock
        if(jpblk(i).ge.0) Cycle
        iblk(itype) = i; jpblk(i) = 0
        nblk(itype) = nblk(itype) + 1
        kblk(itype,nblk(itype)) = i
        m = 1; Exit
       End do
       if(m.ne.0) Return
      end if

! ... everything is full - it is time to generate matrix;
! ... we do it for type with biggest occupation

      jtype = maxloc(nblk,DIM=1)
      Call Collect_coef(jtype)
      if(jtype.eq.itype) Return

! ... find again new block:

      m = 0
      Do i=1,nblock
       if(jpblk(i).ge.0) Cycle
       iblk(itype) = i; jpblk(i) = 0
       nblk(itype) = nblk(itype) + 1
       kblk(itype,nblk(itype)) = i
       m = 1; Exit
      End do
      if(m.eq.0) Stop ' Add_coef_cblock: problems with new block '

      End Subroutine Add_coef_cblock



!======================================================================
      Subroutine Add_c_block(ip,jp,C,j1,j2,j3,j4)
!======================================================================
!     add new data to the sub-list (ip:jp) in module cmdata
!----------------------------------------------------------------------
      Use c_blocks

      Implicit none
      Integer:: ip,jp, j1,j2,j3,j4, k,l,m, i,ii
      Real(8), intent(in) :: C

! ... case of first element in the list:

      if(jp.lt.ip) then
       cdata(ip)=C; k1(ip)=j1; k2(ip)=j2; k3(ip)=j3; k4(ip)=j4
       jp=ip; Return
      end if

! ... search position (k) for new integral

      k=ip; l=jp;
    1 if(k.gt.l) go to 2
      m=(k+l)/2
      if(m.lt.ip.or.m.gt.jp) Stop 'Promlems with m in Add_c_block'
      if    (j1.lt.K1(m)) then;       l = m - 1
      elseif(j1.gt.K1(m)) then;       k = m + 1
      else
       if    (j2.lt.K2(m)) then;      l = m - 1
       elseif(j2.gt.K2(m)) then;      k = m + 1
       else
        if    (j3.lt.K3(m)) then;     l = m - 1
        elseif(j3.gt.K3(m)) then;     k = m + 1
        else
         if    (j4.lt.K4(m)) then;    l = m - 1
         elseif(j4.gt.K4(m)) then;    k = m + 1
         else
          cdata(m)=cdata(m)+C; Return ! the same integral
         end if
        end if
       end if
      end if
      go to 1
    2 Continue

! ... shift the rest of data up:

      if(k.lt.ip.or.k.gt.jp+1) Stop 'Promlems with k in Add_c_block'
      if(jp-ip+1.ge.mblock) Stop 'Promlems with jp in Add_c_block'
      if(jp.ge.k) then
      Do i=jp,k,-1
       ii = i + 1
       cdata(ii)=cdata(i)
       K1(ii)=K1(i); K2(ii)=K2(i); K3(ii)=K3(i); K4(ii)=K4(i)
      End do
      end if

! ... add new integral:

      Cdata(k)=C; K1(k)=j1; K2(k)=j2; K3(k)=j3; K4(k)=j4; jp=jp+1

      End Subroutine Add_c_block

!======================================================================
      Subroutine Merge_c_block(nn, ip, jp, nc)
!======================================================================
!     merge the different blocks of data in MODULE 'cmdata'
!     nn    - number of blocks
!     ip(.) - pointer for begining of block .
!     jp(.) - pointer for end of block .
!     nc    - number of result coeff's
!     all coefficients < EPS_cdata are ignored
!----------------------------------------------------------------------
      Use c_blocks

      Implicit none
      Integer :: nn,nc
      Integer :: ip(*),jp(*)
      Integer :: i,ii, j,jj, m,mm

      nc = 0

! ... choose the non-empty block

       mm=0
       Do m=1,nn; if(JP(m).le.0) Cycle; mm=m; Exit;  End do
       if(mm.eq.0) Return

! ...  main loop ...

    1 Continue

! ...  compare integrals in different blocks and merge the coefficients
! ...  in case of equal integrals (nn > 1)

       Do ii=1,nn-1;  if(JP(ii).le.0) Cycle; i=IP(ii)
        Do jj=ii+1,nn; if(JP(jj).le.0) Cycle; j=IP(jj)
         if(K1(i).ne.K1(j)) Cycle
         if(K2(i).ne.K2(j)) Cycle
         if(K3(i).ne.K3(j)) Cycle
         if(K4(i).ne.K4(j)) Cycle
         CDATA(i) = CDATA(i) + CDATA(j); CDATA(j) = 0.d0
         mm=jj; go to 2
        End do
       End do

! ...  choose the minimum K1, then K2, then K3, then K4

       j=IP(mm)
       Do m=1,nn; if(JP(m).eq.0) Cycle; i=IP(m)
        if    (K1(i).lt.K1(j)) then;   mm=m; j=i
        elseif(K1(i).gt.K1(j)) then;   Cycle
        else
         if    (K2(i).lt.K2(j)) then;  mm=m; j=i
         elseif(K2(i).gt.K2(j)) then;  Cycle
         else
          if    (K3(i).lt.K3(j)) then; mm=m; j=i
          elseif(K3(i).gt.K3(j)) then; Cycle
          elseif(K4(i).lt.K4(j)) then; mm=m; j=i
          end if
         end if
        end if
       End do

! ...  mark the chosen coefficient

       if(abs(CDATA(j)).gt.eps_cdata) then; nc=nc+1; IPT(nc)=IP(mm); end if

! ...  choose next data

    2  IP(mm) = IP(mm) + 1
       if(IP(mm).le.JP(mm)) then
        go to 1
       else
        JP(mm)=0
        Do m=1,nn; if(JP(m).gt.0) then; mm=m; go to 1; end if; End do
       end if

      End Subroutine Merge_c_block


!======================================================================
      Subroutine Collect_coef(itype)
!======================================================================
!     merge data and generate the interaction matrix for given type
!     External call: gen_matrix
!----------------------------------------------------------------------
      Use c_blocks, nc => ncdata

      Implicit none
      Integer, intent(in) :: itype
      Integer :: i,j,k,n

! ... prepare the data:

      n = nblk(itype)                    ! number of blocks
      i = kblk(itype,1)                  ! first block
      if(n.eq.1.and.jpblk(i).le.0) n=0

      if(n.le.0) then                    ! nothing to do
       nc=0; Return

      elseif(n.eq.1) then                ! simple case - one block
       nc = jpblk(i)-ipblk(i) + 1
       j = ipblk(i)-1
       Do i = 1,nc;  IPT(i)=j+i;  End do

      else                               ! need merging data from
                                         ! different blocks
       Do i = 1,n
        j = kblk(itype,i)
        ipi(i) = ipblk(j)
        ipj(i) = jpblk(j)
       End do
       Call Merge_c_block(n, ipi, ipj, nc)

       ! .. release the blocks:

       Do i=1,n; jpblk(kblk(itype,i))=-1;  End do

      end if

! ... re-assign the first block for given itype:

      i=kblk(itype,1); iblk(itype)=i; nblk(itype)=1; jpblk(i)=0

! ... record or use the generated coefficiebts:

      Call Do_coef(itype)

      End Subroutine Collect_coef


!======================================================================
      Real(8) Function C_blocks_occupation()      Result(S)
!======================================================================
      Use c_blocks

      Implicit none
      Integer :: i

      S = 0.d0
      Do i=1,nblock
       if(jpblk(i).le.0) Cycle
       S = S + jpblk(i)-ipblk(i)
      End do

      S = S / (mblock*nblock)

      End Function C_blocks_occupation





