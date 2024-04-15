!======================================================================
      Subroutine plot_bsw
!======================================================================
! ... Computes and tabular the radial orbitals in all gausian points,
! ... plus border values, for further plots
! ... at small r,  P = P0 * r^gamma (1 + r ...)
!----------------------------------------------------------------------
      Use zconst,      only: c_au
      Use DBS_nuclear, only: nuclear
      Use DBS_grid
      Use DBS_gauss
      Use df_orbitals
      Use dbsr_hf

      Implicit none
      Real(8) :: yp(nv*ks+2,nbf),yq(nv*ks+2,nbf),r(nv*ks+2)
      Integer :: i,j,io,i1,i2,m
      Character(6) :: pel(nbf), qel(nbf)
      Real(8) :: P0, gamma

      if(out_w.eq.0.and.out_plot.eq.0) Return
      AF_plt = trim(name)//BF_plt
      Call Read_apar(inp,'plot',AF_plt)
      Call Read_aarg('plot',AF_plt)
      open(nup,file=AF_plt)

! ... radial points for output:

      m=1; r(1)=t(1)
      Do i=1,nv; Do j=1,ks; m=m+1; R(m)=gr(i,j); End do; End do
      m=m+1; R(m)=t(ns+1)

      yp = 0.d0
      yq = 0.d0
      Do io=1,nbf

       Call Bvalue_bm(ksp,p(1,1,io),yp(1,io),pbsp)
       Call Bvalue_bm(ksq,p(1,2,io),yq(1,io),qbsp)
       pel(io)='p'//ebs(io); Call Clean_a(pel(io))
       qel(io)='q'//ebs(io); Call Clean_a(qel(io))

       gamma = lbs(io) + 1
       if(nuclear.eq.'point') gamma = sqrt (kbs(io)**2 - (z/c_au)**2)
       P0 = 0.d0
       Do i = 1,ksp; j = i + 1
        P0 = P0 + yp(j,io)/r(j)**gamma
       End do
       P0 = P0 / ksp; yp(1,io) = P0
      ! or simply ???
       yp(1,io) = yp(2,io)/r(2)**gamma
      End do

      if(out_plot.gt.0) then
      rewind(nup)
      i2=nbf; i1=max(1,i2-nit+1)
      write(nup,'(100(6x,a6))') 'r',(pel(i),qel(i),i=i1,i2)
      Do i=1,m
       write(nup,'(100(1PE12.4))') r(i),(yp(i,io),yq(i,io),io=i1,i2)
      End do
      end if

      End Subroutine plot_bsw


!======================================================================
      Subroutine GRASP_wfn
!======================================================================
! ... output in GRASP format:
!----------------------------------------------------------------------
      Use zconst,      only: c_au
      Use DBS_nuclear, only: nuclear
      Use DBS_grid
      Use DBS_gauss
      Use df_orbitals
      Use dbsr_hf

      Implicit none

      Integer, parameter :: ng = 540  ! max. number of points in GRASP
      Real(8) :: yp(ng),yq(ng),r(ng)
      Real(8) :: P0, gamma, r_max, RNT,HNT
      Integer :: i,io,m,np,nr
      Real(8), external :: bvalu2

      if(out_w.eq.0) Return

! ... radial points for output:

      if(nuclear.eq.'point') then
       RNT = EXP (-65.0d0/16.0d0) / z
       HNT = 0.5d0**4
       np  = ng
      else
       RNT = 2.d-6
       HNT = 5.d-2
       np  = min(ng,220)
      end if

      r_max = 0.d0
      Do io = 1,nwf
       if(r_max.lt.t(mbs(io)+ks)) r_max=t(mbs(io)+ks)
      End do

      Do
       Do i=1,np
        r(i)=RNT*(exp((i-1)*HNT)-1.d0)
       End do
       if(r(np).gt.r_max) Exit
       HNT = HNT*1.2
      End do
      r(1) = 0.d0

! ... file:

      AF_rwf = trim(name)//BF_rwf
      Call Read_ipar(inp,'w',AF_rwf)
      Call Read_iarg('w',AF_rwf)
      Open(nuw,file=AF_rwf,form='UNFORMATTED')
      rewind(nuw)
      write(nuw) 'G92RWF'

! ... cycle over orbitals

      Do io = 1,nbf; m = mbs(io)

       r_max = t(m+ks)
       yp = 0.d0; if(m.lt.ns) p(m+1:ns,1,io)=0.d0
       yq = 0.d0; if(m.lt.ns) p(m+1:ns,2,io)=0.d0
       Do i = 2,np-1
        yp(i) =  bvalu2 (tp, p(1,1,io), nsp, ksp, r(i), 0)
        yq(i) =  bvalu2 (tq, p(1,2,io), nsq, ksq, r(i), 0)
        nr = i; if(r(i+1).gt.r_max) Exit
       End do
       nr = nr + 1

       gamma = lbs(io) + 1
       if(nuclear.eq.'point') gamma = sqrt (kbs(io)**2 - (z/c_au)**2)
       P0 = yp(2)/r(2)**gamma

       write(nuw) nbs(io),kbs(io),e(io,io),nr
       write(nuw) P0,yp(1:nr),yq(1:nr)
       write(nuw) r(1:nr)

      End do

      Call Write_isodata

      End Subroutine GRASP_wfn


!======================================================================
      Subroutine Write_isodata
!======================================================================
      Use DBS_nuclear
      Use dbsr_hf, only: Etotal

      Implicit none
      Integer :: nu
      Real(8) :: apar, cpar

      nu = 91
      open(nu,file='isodata')

      WRITE(nu,'(a)') 'Atomic number:'
      WRITE(nu,*) atomic_number

      WRITE(nu,'(a)') 'Mass number (integer) :'
      WRITE(nu,*) NINT(atomic_weight)

       apar = a_fermi / fermi_in_cm * bohr_radius_in_cm
       cpar = c_fermi / fermi_in_cm * bohr_radius_in_cm


      WRITE(nu,'(a)') 'Fermi distribution parameter a:'
      WRITE(nu,*) apar, a_fermi

      WRITE(nu,'(a)') 'Fermi distribution parameter c:'
      WRITE(nu,*) cpar, c_fermi

      WRITE(nu,'(a)') 'Mass of nucleus (in amu):'
      WRITE(nu,*) atomic_weight

      WRITE(nu,'(a)') 'Nuclear spin (I) (in units of h / 2 pi):'
      WRITE(nu,*) 0.d0

      WRITE(nu,'(a)') 'Nuclear dipole moment (in nuclear magnetons):'
      WRITE(nu,*) 0.d0

      WRITE(nu,'(a)') 'Nuclear quadrupole moment (in barns):'
      WRITE(nu,*) 0.d0

      write(nu,*)
      write(nu,'(a,E15.5,a)') 'RRMS  =',RRMS,   '  !  root mean squared'
      write(nu,*)
      write(nu,'(a,E15.5,a)') 'TPARM =',t_fermi,'  !  skin thickness'
      write(nu,*)
      write(nu,'(a,i5,10x,a)')'nelc  =',NINT(atomic_number),   '  !  nuumber of electrons in neutral atom'
      write(nu,*)
      write(nu,'(a,E15.5,a)') 'Ebind =',Etotal, '  !  gound state energy in Hartrees'


      End Subroutine Write_isodata
