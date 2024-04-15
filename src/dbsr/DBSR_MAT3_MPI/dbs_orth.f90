!=====================================================================
      Subroutine DBS_ORTH
!=====================================================================
!     Imposing orthogonal conditions for scattering orbitals
!     by updating the interaction matrix according the procedure
!     suggested by M.Bentley J.Phys.B27 (1994) 637.
!
!     Let i-th channel is supposed to be orthogonal to orbital p.
!     It can be achieved by following modification of Hamiltonian:
!
!                H -->  (1 - Bcc') H (1 - cc'B)
!
!                H -->  (1 - a) H (1 - b) = H - aH - Hb + aHB
!
!     aH -> multiplication of row by a :  a Hij, j=1,n
!     Hb -> multiplication of column by b: Hij b, i=1,n
!     where c - full solution vector with all zero elements except
!     the channel i is replaced on B-spline expansion of orbital p
!     primes denote the transportasion.
!
!     It leads to transformation of separate blocks  H(i,j):
!
!     H(i,j) with j <> i -->  (1 - Bpp') H(i,j)
!     H(j,i) with j <> i -->   H(i,j) (1 - pp'B)
!     H(i,i) with j  = i -->  (1 - Bpp') H(i,j) (1 - pp'B)
!     j also includes the block for interaction with perturbers
!---------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Real(8) :: S, vv(ms),ww(ms), aa(ms,ms), bb(ms,ms),cc(ms,ms), x(ms,ms)
      Integer :: i,j,ich,jch,ii,jj,ij,ic,jc, nort
      Integer, external :: IBORT

      Do ich = 1,nch; i = ipch(ich)

       nort = 0; aa = 0.d0;  bb = 0.d0

       Do j = 1,nbf
        if(ipbs(j).gt.0) Cycle
        if(kbs(j).ne.kbs(i)) Cycle
        if(IBORT(i,j).ne.0) Cycle
        nort = nort + 1
        Call Get_pv(j,vv,ns)
        Call Get_qv(j,ww,ns)
        Do ii=1,ms; Do jj=1,ms
         S = ww(ii)*vv(jj); aa(ii,jj) = aa(ii,jj)+S; bb(jj,ii) = bb(jj,ii)+S
        End do; End do
       End do

       if(nort.eq.0) Cycle

       Do jch = 1,nch

        if(idiag.eq.+1.and.ich.ne.jch) Cycle
        if(idiag.eq.-1.and.ich.eq.jch) Cycle

        ij = icc(ich,jch); if(ij.eq.0) Cycle

        if(ich.gt.jch) then
         x = hch(:,:,ij)
         cc = matmul(aa,x);   x = x - cc
        elseif(ich.eq.jch) then
         x = hch(:,:,ij)
         cc = matmul(aa,x);   x = x - cc
         cc = matmul(x,bb);   x = x - cc
        else
         x = hch(:,:,ij)
         cc = matmul(x,bb);   x = x - cc
        end if

        hch(:,:,ij) = x

       End do ! over channels (jch)

       if(npert.gt.0.and.idiag.eq.0) then
        Do ic = 1,npert;  jc = icb(ich,ic); if(jc.eq.0) Cycle
         vv = matmul(hcp(:,jc),bb); hcp(:,jc)=hcp(:,jc)-vv
        End do
       end if

      End do ! over channels (ich)

      End Subroutine DBS_ORTH


!=====================================================================
      Subroutine Pri_orth
!=====================================================================
!     printing the imposed orthogonal conditions
!---------------------------------------------------------------------
      Use dbsr_mat

      Implicit none
      Character(200) :: line
      Integer :: i,j,ich,iline
      Integer, external :: IBORT

      write(pri,'(/a/)') 'Orthogonal conditions:'

      Do ich = 1,nch; i = ipch(ich)

       write(line,'(a5,a4)') ebs(i),' -> '; iline=9

      Do j = 1,nbf
       if(ipbs(j).gt.0) Cycle
       if(kbs(j).ne.kbs(i)) Cycle
       if(IBORT(i,j).ne.0) Cycle

       write(line(iline:),'(a5)') ebs(j); iline=iline+6

       if(iline.ge.194) then
        write(pri,'(a)') trim(line)
        iline=9; line = ' '
       end if
      End do ! over bound orbitals (j)

      if(len_trim(line).gt.9) write(pri,'(a)') trim(line)

      End do ! over channels (ich)

      End Subroutine Pri_orth



