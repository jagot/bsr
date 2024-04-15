!======================================================================
      Subroutine Solv_mat
!======================================================================
!     solves the generalized eigenvalue problem  A x = E C x
!----------------------------------------------------------------------
      Use dbsr_pol
      Use DBS_grid
      Use DBS_orbitals_pq
      Use channel_jj, only: nch,npert,ipch,ELC,ipar,jpar

      Implicit none
      Integer :: i,j,ich,ii,jj,m,k,i1,i2,j1,j2,jb, info
      Integer, external :: Ipointer, IORT
      Real(8), allocatable :: sol(:), aa(:),  cc(:,:), c(:)
      Real(8) :: S,EP, dmat, fvalue, alpha
      Real(8), external :: quadr
      Character(64) :: LAB

      if(allocated(sol)) Deallocate(sol); Allocate(sol(mhm))
      if(allocated(aa) ) Deallocate(aa ); Allocate(aa (mhm))

! ... add orthogonal constraints:

       jj = nhm
!      Do ich=1,nch; i=ipch(ich)
!       Do jb =1,nbf
!        if(kbs(i).ne.kbs(jb)) Cycle
!        if(ipbs(jb).ne.0) Cycle
!        if(IORT(i,jb).ne.0) Cycle
!        jj=jj+1
!        j1 = ipsol(ich-1)+1; j2=ipsol(ich)
!        Do j=j1,j2
!         S = quadr(pq(1,1,jb),bb(1,j),0)
!         hm(jj,j) = S
!         hm(j,jj) = S
!        End do
!       End do
!      End do

! ... additional constraints:

      rewind(nuq)
      Do i=1,nortb
       read(nuq) aa(1:nhm)
       jj=jj+1
       hm(1:nhm,jj) = aa(1:nhm)
       hm(jj,1:nhm) = aa(1:nhm)
      End do

! ... energy part:

      hm(1:nhm,1:nhm) = hm(1:nhm,1:nhm) - E1*om(1:nhm,1:nhm)
      sol = d

! ... solve the equation:

      Call LAP_DGESV(mhm,mhm,1,hm,sol,info)
      if(info.ne.0) Stop 'DBSR_POL: solution failed'

! ... normalize the solution:

      S = zero
      Do i = 1,nhm; Do j = 1,nhm
       S = S + sol(i)*om(i,j)*sol(j)
      End do; End do
      S = sqrt(S);
      sol = sol / S

      write(pri,'(/a/)') 'Solution: '

      write(pri,'(a,f16.8,a)') 'norma  = ',S,'  - normalization constant'

! ... restore the solutions in original B-spline basis:

      v = 0.d0
      Do i = 1,nch; i1=(i-1)*ms+1; i2=i*ms
       j1 = ipsol(i-1)+1; j2=ipsol(i)
       Do j=j1,j2
        v(i1:i2) = v(i1:i2) + sol(j)*bb(1:ms,j)
       End do
      End do
      if(npert.gt.0) v(nch*ms+1:khm)=sol(nsol+1:nhm)

! ... energy:

      rewind(nua)
      Do j=1,nhm; read(nua) hm(1:nhm,j); End do

      EP = 0.0
      Do i = 1,nhm; Do j = 1,nhm
       EP = EP + sol(i)*hm(i,j)*sol(j)
      End do; End do
      Deallocate(hm)

      write(*,*) 'EP = ',EP
      write(pri,'(a,f16.8,a)') 'EP     = ',EP,'  - pseudo-state energy, a.u.'
      write(pri,'(a,f16.8,a)') 'EP-E1  = ',(EP-E1)*27.2113, &
                                               '  - pseudo-state energy, eV'
! ... alpha:

      jot1 = jot1+1
      dmat   = Sum(sol(1:nhm)*d(1:nhm))
      fvalue = 2.d0*dmat*dmat/((kpol+kpol+1)*jot1)*(EP-E1)
      alpha  = 2.d0*dmat*dmat/((EP-E1)*(kpol+kpol+1)*jot1)

      write(*,*) 'alpha = ',alpha,'  kpol =',kpol

      write(pri,*)
      write(pri,'(a,f16.8,a)') 'dmat   = ',dmat,  '  - multiole m.e.'
      write(pri,'(a,f16.8,a)') 'fvalue = ',fvalue,'  - oscillator strength'
      write(pri,'(a,f16.8,a)') 'alpha  = ',alpha, '  - polarizability, a.u.'

! ... check the orthogonality:

      write(pri,'(/a/)') 'orthogonality check:'
      Do ich=1,nch; i=ipch(ich); ii=(ich-1)*ms+1
       Do j = 1,nbf;  if(kbs(i).ne.kbs(j)) Cycle
        if(IORT(i,j).ne.0) Cycle
        S = quadr(pq(1,1,j),v(ii),0)
        write(pri,'(2a6,f13.8)') ebs(i),ebs(j),S
       End do
      End do

      Call Check_nortb(v)

! ... output of solution:

      ii = INDEX(AF_pol,'.'); AF = AF_pol(1:ii)//ALSP
      Open(nur,file=AF,form='UNFORMATTED')

      write(nur) khm, nch, npert, ns, jpar, ipar, 1
      LAB = ELC(1)
      write(nur) 1,LAB
      write(nur) EP,(EP-E1)*27.2113, 1, 1
      write(nur) v(1:khm),sol(1:nhm)

      Deallocate(sol, aa)

      End Subroutine Solv_mat

