!======================================================================
      Real(8) Function DBS_core_pq (nclosed,mbreit)
!======================================================================
!     computes the energy of the common closed shells (called as core);
!     core is supposed to have "nclosed" shells in exact correspondence
!     with first "nclosed" orbitals in module DBS_orbitals_pq.
!     mbreit - flag to include the Breit interaction
!     Program first collects all integrals with coefficients and then
!     evaluate them to minimized call to integrals subroutines
!----------------------------------------------------------------------
      Use DBS_orbitals_pq
      Use rk4_data

      Implicit none
      Integer, intent(in) :: nclosed, mbreit
      Integer :: i,j,k,v, ka,ja,la, kb,jb,lb, int, i1,i2,j1,j2, ii,jj,ij,ji
      Real(8) :: C, ca,cb, CK, S(8)
      Integer, external :: l_kappa,j_kappa
      Real(8), external :: Cjkj, Vp_dhl, SMU, zint
      Real(8), parameter :: zero=0.d0, one=1.d0, two=2.d0, half=0.5d0, &
                            eps_c = 1.d-6
      Integer, parameter :: ibr=2**16

      DBS_core_pq = zero; if(nclosed.le.0) Return

      Call alloc_rk_data(irk); nrk=0

      Do i=1,nclosed;  ii = i*ibr+i
       ka=kbs(i); ja=j_kappa(ka); la=l_kappa(ka); ca=ja+1

       DBS_core_pq = DBS_core_pq + ca*Vp_dhl(ka,pq(1,1,i),ka,pq(1,1,i))

       C=ca*ja/two;  k=0;  Call Add_rk4_data(3,k,ii,ii,C)
                           Call Add_rk4_data(4,k,ii,ii,C)
                           Call Add_rk4_data(5,k,ii,ii,C)
                           Call Add_rk4_data(6,k,ii,ii,C)
       Do k = 2,2*la,2
        C = -Cjkj(ja,k,ja)**2 / two
        Call Add_rk4_data(3,k,ii,ii,C)
        Call Add_rk4_data(4,k,ii,ii,C)
        Call Add_rk4_data(5,k,ii,ii,C)
        Call Add_rk4_data(6,k,ii,ii,C)
       End do

       Do j=i+1,nclosed; jj = j*ibr+j; ij=i*ibr+j; ji=j*ibr+i
        kb=kbs(j); jb=j_kappa(kb); lb=l_kappa(kb); cb=jb+1

        C=ca*cb; k=0; Call Add_rk4_data(3,k,ii,jj,C)
                      Call Add_rk4_data(4,k,ii,jj,C)
                      Call Add_rk4_data(5,k,ii,jj,C)
                      Call Add_rk4_data(6,k,ii,jj,C)

        Do k = iabs(la-lb),la+lb,2
         C = -Cjkj(ja,k,jb)**2
         Call Add_rk4_data(3,k,ij,ji,C)
         Call Add_rk4_data(4,k,ij,ji,C)
         Call Add_rk4_data(5,k,ij,ji,C)
         Call Add_rk4_data(6,k,ij,ji,C)
        End do

      End do; End do

! ... evaluate the integrals:

      C = zero
      Do i = 1,nrk; int=kr1(i); k=kr2(i)
       if(abs(crk(i)).lt.eps_C) Cycle
       i1=kr3(i)/ibr; i2=mod(kr3(i),ibr)
       j1=kr4(i)/ibr; j2=mod(kr4(i),ibr)
       CK = zint(int,i1,j1,i2,j2,k)
       C = C + crk(i)*CK
      End do

      DBS_core_pq = DBS_core_pq + C

!----------------------------------------------------------------------
! ... Breit conttribution:

      if(mbreit.eq.0) Return

      nrk=0; int=7

      Do i=1,nclosed;  ka=kbs(i); ja=j_kappa(ka)
      Do j=i,nclosed;  kb=kbs(j); jb=j_kappa(kb)

       ij=i*ibr+j; ji=j*ibr+i

       Do k = iabs(ja-jb)/2,(ja+jb)/2

        C = -Cjkj(ja,k,jb)**2;   if(i.eq.j) C = C / two
        if(C.eq.zero) Cycle

        Do v = k-1,k+1
         if(SMU(ka,kb,kb,ka,k,v,S).eq.0.d0) Cycle

         Call Add_rk4_data(int,v,ij,ji,C*S(1))
         Call Add_rk4_data(int,v,ji,ij,C*S(2))
         Call Add_rk4_data(int,v,ji,ij,C*S(3))
         Call Add_rk4_data(int,v,ij,ji,C*S(4))
         Call Add_rk4_data(int,v,ij,ij,C*S(5))
         Call Add_rk4_data(int,v,ij,ij,C*S(6))
         Call Add_rk4_data(int,v,ji,ji,C*S(7))
         Call Add_rk4_data(int,v,ji,ji,C*S(8))
        End do

       End do

      End do; End do

! ... evaluate the integrals:

      C = zero
      Do i = 1,nrk; int=kr1(i); k=kr2(i)
       if(abs(crk(i)).lt.eps_C) Cycle
       i1=kr3(i)/ibr; i2=mod(kr3(i),ibr)
       j1=kr4(i)/ibr; j2=mod(kr4(i),ibr)
       if(int.ne.7) Stop 'DBS_core: SK = ?'
       CK = zint(int,i1,j1,i2,j2,k)
       C = C + crk(i)*CK
      End do

      DBS_core_pq = DBS_core_pq + C

      End Function DBS_core_pq
