!======================================================================
!     PROGRAM       H _ T X T
!
!               C O P Y R I G H T -- 2024
!
!     Written by:   Stefanos Carlström
!                   email: stefanos.carlstrom@gmail.com
!
!======================================================================
!  Translate h.nnn -> h.nnn.txt
!======================================================================
!
!  INPUT PARAMETERS:
!
!     klsp1     -  Starting partial wave [default 1]
!     klsp2     -  Final partial wave [default 999]
!
!  INPUT FILES:
!
!     h.nnn     -  Surface data, resolved on partial wave
!
!  OUTPUT FILES:
!
!     h.nnn.txt -  Textual representation of surface data
!
!======================================================================
program h_txt
  Implicit none
  Integer :: klsp1=1, klsp2=999, klsp
  character(256) :: name=""

  Call read_name(name)
  if(name == "?") call print_help

  Call Read_iarg('klsp1', klsp1)
  Call Read_iarg('klsp2', klsp2)

  do klsp = klsp1,klsp2
     call translate(klsp)
  end do
contains
  subroutine print_help
    write(*,*) "Usage: h_txt klsp1= klsp2="
    stop
  end subroutine print_help

  subroutine translate(klsp)
    use file_helper

    Implicit none

    Integer, intent(in) :: klsp
    character(256) :: infilename, outfilename

    type(file) :: in, out
    Integer :: i, j, k

    Integer, External :: Icheck_file

    Integer :: nelc, nz, lrange, km, ntarg, jpar, spar, par, nch, khm, more
    Real(8) :: RA, RB

    Real(8), allocatable :: etarg(:)
    Integer, allocatable :: jtarg(:), starg(:), ptarg(:), nconat(:), iptar(:)

    Integer, allocatable :: lch(:), kch(:)
    Real(8), allocatable :: cf(:,:,:), eval(:), wmat(:,:)

    Character(len=*), parameter :: intfmt = '(a10," = ",i5)', &
         realfmt = '(a10," = ",e26.16e3)'

    write(infilename, '("h.",i0.3)') klsp
    write(outfilename, '(a,".txt")') trim(infilename)
    if(Icheck_file(trim(infilename)) == 0) return

    write(*,'("Translating partial-wave surface data from ",a," to ",a)') &
         trim(infilename), trim(outfilename)

    open(newunit=in%unit, file=trim(infilename), form='UNFORMATTED')
    read(in%unit) nelc, nz, lrange, km, ntarg, RA, RB
    allocate(etarg(ntarg),jtarg(ntarg),starg(ntarg),ptarg(ntarg),nconat(ntarg))

    read(in%unit) etarg
    read(in%unit) jtarg
    read(in%unit) starg, ptarg

    ! Skip the (non-existent) Buttle correction
    read(in%unit)

    read(in%unit) jpar, spar, par, nch, khm, more

    allocate(lch(nch),kch(nch),iptar(nch), &
         cf(nch,nch,km),eval(khm),wmat(nch,khm))

    read(in%unit) nconat
    call channel_target_mapping(nconat, iptar)
    read(in%unit) lch, kch
    read(in%unit) cf
    read(in%unit) eval
    read(in%unit) wmat

    open(newunit=out%unit, file=trim(outfilename), form='FORMATTED')
    write(out%unit,intfmt) "nelc", nelc
    write(out%unit,intfmt) "nz", nz
    write(out%unit,intfmt) "lrange", lrange
    write(out%unit,intfmt) "km", km
    write(out%unit,intfmt) "ntarg", ntarg
    write(out%unit,realfmt) "RA", RA
    write(out%unit,realfmt) "RB", RB

    write(out%unit,intfmt) "jpar", jpar
    write(out%unit,intfmt) "spar", spar
    write(out%unit,intfmt) "par", par
    write(out%unit,intfmt) "nch", nch
    write(out%unit,intfmt) "khm", khm

    write(out%unit,*)
    write(out%unit,'("# Target states:")')
    write(out%unit,'("#  Q = L for LS-coupled calculations")')
    write(out%unit,'("#  Q = 2J for LSJ- and JJ-coupled calculations")')
    write(out%unit,'("#",a5,a26,4a5)') "i", "E [Ha]", "Q", "S", "par", "nch"
    do i=1,ntarg
       write(out%unit,'(" ",i5,e26.16e3,4i5)') i, &
            etarg(i), jtarg(i), starg(i), ptarg(i), nconat(i)
    end do

    write(out%unit,*)
    write(out%unit,'("# Channels:")')
    write(out%unit,'("#",5a6)') "i", "l", "kappa", "2j", "targ"
    do i=1,nch
       write(out%unit,'(" ",5i6)') i, lch(i), kch(i), 2*abs(kch(i))-1, iptar(i)
    end do

    write(out%unit,*)
    write(out%unit,'("# Channel asymptotic multipoles:")')
    write(out%unit,'("#",3a5,a26)') "i", "j", "k", "cf(i,j,k)"
    do i=1,nch
       do j=1,nch
          do k=1,km
             write(out%unit,'(" ",3i5,e26.16e3)') i, j, k, cf(i,j,k)
          end do
       end do
    end do

    write(out%unit,*)
    write(out%unit,'("# R-matrix poles and surface amplitudes:")')
    write(out%unit, '("#",a5,a26)', advance='no') "i", "E [Ha]"
    do j=1,nch
       write(out%unit, '(i26)', advance='no') j
    end do
    write(out%unit,*)
    do i=1,khm
       write(out%unit, '(" ",i5,e26.16e3)', advance='no') i, eval(i)
       do j=1,nch
          write(out%unit, '(e26.16e3)', advance='no') wmat(j,i)
       end do
       write(out%unit,*)
    end do

    deallocate(etarg,jtarg,starg,ptarg,nconat,iptar, &
         lch, kch, cf, eval, wmat)
  end subroutine translate

  subroutine channel_target_mapping(nconat, iptar)
    implicit none

    Integer, intent(in) :: nconat(:)
    Integer, intent(out) :: iptar(:)

    Integer :: it, i1, i2, i, ntarg

    ntarg = size(nconat, 1)

    ! Taken from sum_hh_jj.f90
    iptar = 0
    i1=1
    do it=1,ntarg
       if(nconat(it).eq.0) cycle
       i2=i1+nconat(it)-1
       do i=i1,i2
          iptar(i)=it
       end do
       i1=i2+1
    end do
  end subroutine channel_target_mapping
end program h_txt
