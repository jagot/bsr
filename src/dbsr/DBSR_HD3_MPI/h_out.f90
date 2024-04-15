!======================================================================
      Subroutine H_OUT
!======================================================================
!     Output results in H.DAT file (unit 'nuh') for further asymptotic
!     calculations
!----------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Integer :: NCONAT(ntarg)
      Integer :: i,j,k,km, LRANG2,NPTY
      Real(8), allocatable :: acf(:,:,:)   ! asymptotic coefficients

! ... read asympotic coefficients from the dbsr_mat.nnn file:

      read(nui) km
      Allocate(acf(nch,nch,0:km))
      read(nui) acf

      RA = tmax; RB = zero

! ... open output h.nnn file:

      i = INDEX(AF_h,'.'); AF = AF_h(1:i)//ALSP
      Open(nuh,file=AF,form='UNFORMATTED')

      LRANG2 = maxval(lch) + 1

! ... basic parameter:

      NCONAT=0
      write(nuh) nelc,nz,LRANG2,km,ntarg,RA,RB
      write(nuh) E_exp(1:ntarg) ! etarg(1:ntarg)
      write(nuh) jtarg(1:ntarg)
      write(nuh) (0,i=1,ntarg)

! ... Buttle corrections - don't Used in BSR !

      write(nuh) (0.d0,0.d0,0.d0,i=1,LRANG2)

! ... partial wave parameters:

      if(ipar.eq.+1) NPTY=0
      if(ipar.eq.-1) NPTY=1

      write(nuh) jpar,0,NPTY,nch,khm,0

! ... number of channels, coupled to each target state:

      Do i=1,nch; j=iptar(i); NCONAT(j)=NCONAT(j)+1; End do
      write(nuh) NCONAT

! ... angular momentum of continuum electron:

      write(nuh) (lch(i),i=1,nch),(kch(i),i=1,nch)

! ... asymptotic coefficients:

      write(nuh) (((acf(i,j,k),i=1,nch),j=1,nch),k=1,km)
      Deallocate(acf)

! ... R_matrix poles in dereasing order:

      write(nuh) (eval(i),i=khm,1,-1)

! ... surface amplitudes:

      write(nuh) ((WMAT(i,j),i=1,nch),j=khm,1,-1)

      Close(nuh)

      End Subroutine H_OUT


!======================================================================
      Subroutine H_OUT_exp
!======================================================================
!     Output results in H.DAT file (unit 'nuh') for further asymptotic
!     calculations
!----------------------------------------------------------------------
      Use dbsr_hd

      Implicit none
      Integer :: NCONAT(ntarg), npch(nch)
      Integer :: i,j,k,km,it,ich, LRANG2,NPTY
      Real(8), allocatable :: acf(:,:,:)   ! asymptotic coefficients

! ... read asympotic coefficients from the dbsr_mat.nnn file:

      read(nui) km
      Allocate(acf(nch,nch,0:km))
      read(nui) acf

      RA = tmax; RB = zero

! ... open output h.nnn file:

      i = INDEX(AF_h,'.'); AF = AF_h(1:i)//ALSP
      Open(nuh,file=AF,form='UNFORMATTED')

      LRANG2 = maxval(lch) + 1

! ... basic parameter:

      NCONAT=0
      write(nuh) nelc,nz,LRANG2,km,ntarg,RA,RB
      write(nuh) (E_exp(ip_exp(i)),i=1,ntarg)
      write(nuh) (jtarg(ip_exp(i)),i=1,ntarg)
      write(nuh) (0,i=1,ntarg)

! ... Buttle corrections - don't Used in BSR !

      write(nuh) (0.d0,0.d0,0.d0,i=1,LRANG2)

! ... partial wave parameters:

      if(ipar.eq.+1) NPTY=0
      if(ipar.eq.-1) NPTY=1

      write(nuh) jpar,0,NPTY,nch,khm,0

! ... number of channels, coupled to each target state:

      NCONAT=0
      Do i=1,ntarg; it=ip_exp(i)
       Do ich=1,nch; if(iptar(ich).ne.it) Cycle
        NCONAT(i) = NCONAT(i) + 1
       End do
      End do
      write(nuh) NCONAT

! ... new channel order:

      k = 0
      Do i=1,ntarg; it=ip_exp(i)
       Do ich=1,nch; if(iptar(ich).ne.it) Cycle
        k=k+1; npch(k)=ich
       End do
      End do
      if(k.ne.nch) Stop 'Problems with npch in H_OUT_exp'

! ... angular momentum of continuum electron:

      write(nuh) (lch(npch(i)),i=1,nch),(kch(npch(i)),i=1,nch)

! ... asymptotic coefficients:

      write(nuh) (((acf(npch(i),npch(j),k),i=1,nch),j=1,nch),k=1,km)
      Deallocate(acf)

! ... R_matrix poles in dereasing order:

      write(nuh) (eval(i),i=khm,1,-1)

! ... surface amplitudes:

      write(nuh) ((WMAT(npch(i),j),i=1,nch),j=khm,1,-1)

      Close(nuh)

      End Subroutine H_OUT_exp


