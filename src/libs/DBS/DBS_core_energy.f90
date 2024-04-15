!====================================================================
      Module DBS_core_energy
!====================================================================
!     core-associated  parameters for core energy
!--------------------------------------------------------------------
      Implicit none

      Integer :: nbk = 0       ! current number of coefficients
      Integer :: mbk = 0       ! maximum dimension
      Integer :: ibk = 2**10   ! initial dimension
      Integer :: ibi = 2**15   ! packing basis
      Real(8) :: eps_c = 1.d-10

! ... coefficients and their attributes:

      Real(8), allocatable :: cbk(:)
      Integer, allocatable :: kr1(:),kr2(:),kr3(:),kr4(:)

      End Module DBS_core_energy


!======================================================================
      Subroutine Alloc_core_coefs(m)
!======================================================================
! ... allocate, deallocate, or reallocate the data
!----------------------------------------------------------------------
      Use DBS_core_energy

      Implicit none
      Integer :: m
      Integer, allocatable :: iar(:)
      Real(8), allocatable :: rar(:)

      if(m.le.0) then
       if(allocated(cbk)) Deallocate (cbk,kr1,kr2,kr3,kr4)
       nbk = 0; mbk = 0
      elseif(.not.allocated(cbk)) then
       mbk = m; nbk = 0
       Allocate(cbk(mbk),kr1(mbk),kr2(mbk),kr3(mbk),kr4(mbk))
      elseif(m.le.mbk) then
       Return
      elseif(nbk.eq.0) then
       Deallocate (cbk,kr1,kr2,kr3,kr4)
       mbk = m
       Allocate(cbk(mbk),kr1(mbk),kr2(mbk),kr3(mbk),kr4(mbk))
      else
       Allocate(rar(nbk))
       rar=cbk; Deallocate(cbk); Allocate(cbk(m)); cbk(1:nbk)=rar
       Deallocate(rar)
       Allocate(iar(nbk))
       iar=kr1; Deallocate(kr1); Allocate(kr1(m)); kr1(1:nbk)=iar
       iar=kr2; Deallocate(kr2); Allocate(kr2(m)); kr2(1:nbk)=iar
       iar=kr3; Deallocate(kr3); Allocate(kr3(m)); kr3(1:nbk)=iar
       iar=kr4; Deallocate(kr4); Allocate(kr4(m)); kr4(1:nbk)=iar
       Deallocate(iar)
       mbk = m
      end if

      End Subroutine Alloc_core_coefs


!======================================================================
      Subroutine Add_core_coef(int,kpol,i1,i2,i3,i4,C)
!======================================================================
!     add new data to the list
!----------------------------------------------------------------------
      Use DBS_core_energy

      Implicit none
      Integer, intent(in) ::  int,kpol,i1,i2,i3,i4
      Real(8), intent(in) :: C
      Integer :: i,k,l,m, k1,k2,k3,k4

      if(mbk.eq.0) Call Alloc_core_coefs(ibk)

      k1=int; k2=kpol; k3=i1*ibi+i3; k4=i2*ibi+i4

! ... search position (k) for new integral

      k=1; l=nbk
    1 if(k.gt.l) go to 2
      m=(k+l)/2
      if    (k1.lt.kr1(m)) then;       l = m - 1
      elseif(k1.gt.kr1(m)) then;       k = m + 1
      else
       if    (k2.lt.kr2(m)) then;      l = m - 1
       elseif(k2.gt.kr2(m)) then;      k = m + 1
       else
        if    (k3.lt.kr3(m)) then;     l = m - 1
        elseif(k3.gt.kr3(m)) then;     k = m + 1
        else
         if    (k4.lt.kr4(m)) then;    l = m - 1
         elseif(k4.gt.kr4(m)) then;    k = m + 1
         else
          cbk(m)=cbk(m)+C
          Return ! the same integral
         end if
        end if
       end if
      end if
      go to 1
    2 Continue

! ... shift the rest data up:

      Do i=nbk,k,-1
       m = i + 1
       cbk(m)=cbk(i)
       kr1(m)=kr1(i); kr2(m)=kr2(i); kr3(m)=kr3(i); kr4(m)=kr4(i)
      End do

! ... add new integral:

      cbk(k)=C; kr1(k)=k1; kr2(k)=k2; kr3(k)=k3; kr4(k)=k4; nbk=nbk+1
      if(nbk.eq.mbk) Call Alloc_core_coefs(mbk+ibk)

      End Subroutine Add_core_coef


!======================================================================
      Subroutine Get_core_coef (ncore,mbreit,kbs)
!======================================================================
!     define the angular coefficients
!----------------------------------------------------------------------
      Use DBS_core_energy

      Implicit none
      Integer, intent(in) :: ncore, mbreit, kbs(*)
      Integer :: i,j,k,v, ka,ja,la, kb,jb,lb, int, i1,i2,j1,j2
      Real(8) :: C, ca,cb, CK, S(8)
      Integer, external :: l_kappa,j_kappa
      Real(8), external :: Cjkj, SMU

      Real(8), parameter :: zero=0.d0, one=1.d0, two=2.d0, half=0.5d0

      Do i=1,ncore
       ka=kbs(i); ja=j_kappa(ka); la=l_kappa(ka); ca=ja+1

       Call Add_core_coef(0,0,i,i,i,i,ca)

       C=ca*ja/two;  k=0;  Call Add_core_coef(1,k,i,i,i,i,C)
                           Call Add_core_coef(2,k,i,i,i,i,C)
                           Call Add_core_coef(3,k,i,i,i,i,C)
                           Call Add_core_coef(4,k,i,i,i,i,C)
       Do k = 2,2*la,2
        C = -Cjkj(ja,k,ja)**2 / two
        Call Add_core_coef(1,k,i,i,i,i,C)
        Call Add_core_coef(2,k,i,i,i,i,C)
        Call Add_core_coef(3,k,i,i,i,i,C)
        Call Add_core_coef(4,k,i,i,i,i,C)
       End do

       Do j=i+1,ncore
        kb=kbs(j); jb=j_kappa(kb); lb=l_kappa(kb); cb=jb+1

        C=ca*cb; k=0; Call Add_core_coef(1,k,i,j,i,j,C)
                      Call Add_core_coef(2,k,i,j,i,j,C)
                      Call Add_core_coef(3,k,i,j,i,j,C)
                      Call Add_core_coef(4,k,i,j,i,j,C)

        Do k = iabs(la-lb),la+lb,2
         C = -Cjkj(ja,k,jb)**2
         Call Add_core_coef(1,k,i,j,j,i,C)
         Call Add_core_coef(2,k,i,j,j,i,C)
         Call Add_core_coef(3,k,i,j,j,i,C)
         Call Add_core_coef(4,k,i,j,j,i,C)
        End do

      End do; End do

! ... Breit conttribution:

      if(mbreit.eq.0) Return

      int=5
      Do i=1,ncore;  ka=kbs(i); ja=j_kappa(ka)
      Do j=i,ncore;  kb=kbs(j); jb=j_kappa(kb)

       Do k = iabs(ja-jb)/2,(ja+jb)/2

        C = -Cjkj(ja,k,jb)**2;   if(i.eq.j) C = C / two
        if(C.eq.zero) Cycle

        Do v = k-1,k+1
         if(SMU(ka,kb,kb,ka,k,v,S).eq.0.d0) Cycle
         Call Add_core_coef(int,v,i,j,j,i,C*S(1))
         Call Add_core_coef(int,v,j,i,i,j,C*S(2))
         Call Add_core_coef(int,v,j,i,i,j,C*S(3))
         Call Add_core_coef(int,v,i,j,j,i,C*S(4))
         Call Add_core_coef(int,v,i,i,j,j,C*S(5))
         Call Add_core_coef(int,v,i,i,j,j,C*S(6))
         Call Add_core_coef(int,v,j,j,i,i,C*S(7))
         Call Add_core_coef(int,v,j,j,i,i,C*S(8))
        End do

       End do

      End do; End do

      End Subroutine Get_core_coef


!======================================================================
      Real(8) Function Ecore_dbs (ncore,mbreit,kbs)
!======================================================================
!     compute energy of the common closed shells (core)
!----------------------------------------------------------------------
      Use DBS_core_energy

      Implicit none
      Integer, intent(in) :: ncore, mbreit, kbs(*)
      Integer :: i, int,k, i1,i2,j1,j2
      Real(8), external :: zint_pq

      Ecore_dbs = 0.d0;  if(ncore.eq.0) Return

      if(nbk.eq.0) Call Get_core_coef (ncore,mbreit,kbs)

! ... evaluate the integrals:

      Do i = 1,nbk; int=kr1(i); k=kr2(i)
       if(abs(cbk(i)).lt.eps_C) Cycle
       i1=kr3(i)/ibi; i2=mod(kr3(i),ibi)
       j1=kr4(i)/ibi; j2=mod(kr4(i),ibi)
       Ecore_dbs = Ecore_dbs + cbk(i)*zint_pq(int,i1,j1,i2,j2,k)
      End do
      End Function Ecore_dbs
