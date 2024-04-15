!======================================================================
!     PROGRAM       D B S R _ C I
!
!               C O P Y R I G H T -- 2015
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!     A CONFIGURATION INTERACTION program in the DIRAC-COULOMB-BREIT
!     approach including the case of NON-ORTHOGONAL radial orbitals.
!     It is a part of the DBSR project, where B-spline is used for
!     represantation of radial orbitals.
!----------------------------------------------------------------------
!     INPUT:  command-line arguments
!             name.inp_ci - input parameters (optional)
!             name.c      - list of configurations in the GRASP format
!             name.bsw    - radial functions in B-spline representation
!             name.bnk    - data bank for angular coefficients
!
!     OUTPUT: name.log_ci - running information
!             name.j      - expansions and energies for CI solutions
!             name.d      - debug information (optional)
!----------------------------------------------------------------------
!     The parameters in command line has structure:   param=value
!     except 'name' of case - should be given as the first argument
!----------------------------------------------------------------------
!     List of main parameters with default values:
!     (all arguments except "name" are optional)
!
!     name        - name of case
!     mbreit=1    - include (1) or not (0) Breit oprerators
!     msol=0      - how many eigensolution you need (0 -> all)
!     nzero=0     - zero-order dimension (0 -> all configurations)
!     debug=0     - if /= 0, file name.d contain debug information
!
!     eps_c    = 1.0D-10  -  tolerance for coefficients
!     eps_det  = 1.0D-10  -  tolerance for determinant overlaps
!     eps_ovl  = 1.0D-08  -  tolerance for configuration overlaps
!----------------------------------------------------------------------
!     Example of calling:
!
!     dbsr_ci name
!     dbsr_ci name jmin=1 jmax=3 mbreit=0 msol=1
!----------------------------------------------------------------------
      Use dbsr_ci
      Use conf_jj
      Use DBS_orbitals_pq

      Implicit none
      Integer :: j
      Real(8) :: t1,t2,t3
      Real(8), external :: DBS_core_pq

      Call CPU_time(t1)

! ... read parameters

      Call Read_arg

! ... check input files: name.c, name.bsw, name.bnk

      Call Read_data

! ... define number of J-blocks (for each total momentum J)

      Call Def_jblocks

! ... define the fine-turn parameters for configuration energies

      Call Read_shift

! ... define the semi-empirical corrections for integrals if any

      Call Read_int_corr

! ... define core energy:

      Call CPU_time(t2)
      Ecore = DBS_core_pq(ncore,mbreit)
      Call CPU_time(t3)
      write(pri,'(/a,F16.8)') 'Ecore  =',ECORE
      write(pri,'(/a,T30,f8.2,a)') 'CORE:', (t3-t2)/60, ' min'
      write(  *,'(/a,T30,f8.2,a)') 'CORE:', (t3-t2)/60, ' min'

! ... prepare one-electron integrals:

      Call CPU_time(t2)
      Call Gen_DHL_core(nclosed,mbreit,0)
      Call CPU_time(t3)
      write(pri,'(/a,T30,f8.2,a)') 'L_integrals:', (t3-t2)/60, ' min'
      write(  *,'(/a,T30,f8.2,a)') 'L_integrals:', (t3-t2)/60, ' min'

! ... perform calculation for each J-value

      write (nuj,'(a,i7)') 'ncfg = ',ncfg
      write (nuj,'(/a)') 'Solutions:'

      Do j = 1,njbl; jot = JJc(j)

       ncj = Jncfg(j)
       if(njbl.gt.1.or.nzero.le.0.or.nzero.gt.ncj) nzero=ncj
       meiv = msol
       if(meiv.le.0.or.meiv.gt.nzero) meiv=nzero

       write(*,'(/a,i5 )')  'jot  =',jot
       write(*,'( a,3i5)')  'ncfg =',ncj
       write(*,'( a,i5 )')  'nzero=',nzero
       write(*,'( a,i5 )')  'meiv =',meiv

       if(allocated(HM)) Deallocate(HM,SM,DM,EVAL)
       Allocate(HM(ncj,nzero),SM(ncj,nzero),EVAL(nzero),DM(ncj))

       Call CPU_time(t2)
       Call Setup_matrix(JTc1(j),JTc2(j))
       Call CPU_time(t3)
       write(pri,'(a,T30,f8.2,a)') 'Matrix set-up:', (t3-t2)/60, ' min'
       write(  *,'(a,T30,f8.2,a)') 'Matrix set-up:', (t3-t2)/60, ' min'

       if(check_c.gt.0) then
        if(njbl.gt.1) &
         Stop 'Stop: number of J-blocks > 0 with option check_c > 0'
        Call Check_cfile
        Stop ' '
       end if

       Call CPU_time(t2)
       Call DIAG(JTc1(j)-1)
       Call CPU_time(t3)
       write(pri,'(a,T30,f8.2,a)') 'Diagonalization:', (t3-t2)/60, ' min'
       write(  *,'(a,T30,f8.2,a)') 'Diagonalization:', (t3-t2)/60, ' min'

      End do  ! over J-blocks

! ... total number of solutions recorded in j-file:

      write(nuj,'(a,i7)') 'nsol = ',nsol
      Close(nuj)

! ... timing:

      Call CPU_time(t2)
      write(pri,'(/a,T30,f8.2,a/)') 'Time of calculations =', (t2-t1)/60,' min'
      write(  *,'(/a,T30,f8.2,a/)') 'Time of calculations =', (t2-t1)/60,' min'

      End ! Program dbsr_ci


