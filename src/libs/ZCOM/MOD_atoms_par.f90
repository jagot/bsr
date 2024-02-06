!======================================================================
      Module atoms_par
!======================================================================
!     atomic parameters according to periodic table
!----------------------------------------------------------------------
!     Nuclear radii are taken from
!     I. Angeli, K. P. Marinova, Atomic Data and Nuclear
!     Data Tables {\bf 99}, 69-95 (2013).
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!     If no data for the element are available, we use
!     the empirical formula for Z < 91 from:
!     W.R. Johnson, G.Soff, Atomic Data and Nudear Data
!                   Tables {\bf 33}, p.405 (1985)
!     and for Z> 90 from:
!     I. Goidenko, Private communication.
!----------------------------------------------------------------------
      Implicit none
      Integer, parameter :: n_atoms = 104
      Type atomic
       INTEGER :: an
       CHARACTER(2) :: symbol
       CHARACTER(4) :: core
       CHARACTER(40) :: conf
       REAL(8) :: weight
       REAL(8) :: rrms
      End type atomic

      Type(atomic), dimension(n_atoms), parameter, public :: atoms = (/  &
       atomic(   1,  'H ',  '    ',  '1s(1)',                   1.d0, 0.8783d0), &
       atomic(   2,  'He',  '    ',  '1s(2)',                   4.d0, 1.6755d0), &
       atomic(   3,  'Li',  '[He]',  '2s(1)',                   6.d0, 2.5890d0), &
       atomic(   4,  'Be',  '[He]',  '2s(2)',                   9.d0, 2.5190d0), &
       atomic(   5,  'B ',  '[Be]',  '2p-(1)',                 10.d0, 2.4277d0), &
       atomic(   6,  'C ',  '[Be]',  '2p-(1)2p(1)',            12.d0, 2.4702d0), &
       atomic(   7,  'N ',  '[Be]',  '2p-(1)2p(2)',            14.d0, 2.5582d0), &
       atomic(   8,  'O ',  '[Be]',  '2p-(1)2p(3)',            16.d0, 2.6991d0), &
       atomic(   9,  'F ',  '[Be]',  '2p-(2)2p(3)',            19.d0, 2.8976d0), &
       atomic(  10,  'Ne',  '[Be]',  '2p-(2)2p(4)',            20.d0, 3.0055d0), &
       atomic(  11,  'Na',  '[Ne]',  '3s(1)',                  23.d0, 2.9936d0), &
       atomic(  12,  'Mg',  '[Ne]',  '3s(2)',                  24.d0, 3.0570d0), &
       atomic(  13,  'Al',  '[Mg]',  '3p-(1)',                 27.d0, 3.0610d0), &
       atomic(  14,  'Si',  '[Mg]',  '3p-(1)3p(1)',            28.d0, 3.1224d0), &
       atomic(  15,  'P ',  '[Mg]',  '3p-(1)3p(2)',            31.d0, 3.1889d0), &
       atomic(  16,  'S ',  '[Mg]',  '3p-(1)3p(3)',            32.d0, 3.2611d0), &
       atomic(  17,  'Cl',  '[Mg]',  '3p-(2)3p(3)',            35.d0, 3.3654d0), &
       atomic(  18,  'Ar',  '[Mg]',  '3p-(2)3p(4)',            38.d0, 3.4028d0), &
       atomic(  19,  'K ',  '[Ar]',  '4s(1)',                  39.d0, 3.4349d0), &
       atomic(  20,  'Ca',  '[Ar]',  '4s(2)',                  40.d0, 3.4776d0), &
       atomic(  21,  'Sc',  '[Ar]',  '3d-(1)4s(2)',            45.d0, 3.5459d0), &
       atomic(  22,  'Ti',  '[Ar]',  '3d-(2)4s(2)',            48.d0, 3.5921d0), &
       atomic(  23,  'V ',  '[Ar]',  '3d-(3)4s(2)',            51.d0, 3.6002d0), &
       atomic(  24,  'Cr',  '[Ar]',  '3d-(4)3d(1)4s(1)',       52.d0, 3.6452d0), &
       atomic(  25,  'Mn',  '[Ar]',  '3d-(4)3d(1)4s(2)',       55.d0, 3.7057d0), &
       atomic(  26,  'Fe',  '[Ar]',  '3d-(4)3d(2)4s(2)',       56.d0, 3.7377d0), &
       atomic(  27,  'Co',  '[Ar]',  '3d-(4)3d(3)4s(2)',       59.d0, 3.7875d0), &
       atomic(  28,  'Ni',  '[Ar]',  '3d-(4)3d(4)4s(2)',       60.d0, 3.8118d0), &
       atomic(  29,  'Cu',  '[Ar]',  '3d-(4)3d(6)4s(1)',       63.d0, 3.8823d0), &
       atomic(  30,  'Zn',  '[Ar]',  '3d-(4)3d(6)4s(2)',       66.d0, 3.9491d0), &
       atomic(  31,  'Ga',  '[Zn]',  '4p-(1)',                 69.d0, 3.9973d0), &
       atomic(  32,  'Ge',  '[Zn]',  '4p-(2)',                 74.d0, 4.0742d0), &
       atomic(  33,  'As',  '[Zn]',  '4p-(2)4p(1)',            75.d0, 4.0968d0), &
       atomic(  34,  'Se',  '[Zn]',  '4p-(2)4p(2)',            80.d0, 4.1400d0), &
       atomic(  35,  'Br',  '[Zn]',  '4p-(2)4p(3)',            79.d0, 4.1629d0), &
       atomic(  36,  'Kr',  '[Zn]',  '4p-(2)4p(4)',            86.d0, 4.1835d0), &
       atomic(  37,  'Rb',  '[Kr]',  '5s(1)',                  87.d0, 4.1989d0), &
       atomic(  38,  'Sr',  '[Kr]',  '5s(2)',                  88.d0, 4.2240d0), &
       atomic(  39,  'Y ',  '[Kr]',  '4d-(1)5s(2)',            89.d0, 4.2430d0), &
       atomic(  40,  'Zr',  '[Kr]',  '4d-(2)5s(2)',            90.d0, 4.2694d0), &
       atomic(  41,  'Nb',  '[Kr]',  '4d-(4)5s(1)',            93.d0, 4.3240d0), &
       atomic(  42,  'Mo',  '[Kr]',  '4d-(4)4d(1)5s(1)',       92.d0, 4.3151d0), &
       atomic(  43,  'Tc',  '[Kr]',  '4d-(4)4d(1)5s(2)',       97.d0, 0.0000d0), &  ! ???
       atomic(  44,  'Ru',  '[Kr]',  '4d-(4)4d(3)5s(1)',      104.d0, 4.5098d0), &
       atomic(  45,  'Rh',  '[Kr]',  '4d-(4)4d(4)5s(1)',      103.d0, 4.4945d0), &
       atomic(  46,  'Pd',  '[Kr]',  '4d-(4)4d(6)',           108.d0, 4.5563d0), &
       atomic(  47,  'Ag',  '[Kr]',  '4d-(4)4d(6)5s(1)',      109.d0, 4.5638d0), &
       atomic(  48,  'Cd',  '[Kr]',  '4d-(4)4d(6)5s(2)',      114.d0, 4.6087d0), &
       atomic(  49,  'In',  '[Cd]',  '5p-(1)',                115.d0, 4.6156d0), &
       atomic(  50,  'Sn',  '[Cd]',  '5p-(2)',                120.d0, 4.6519d0), &
       atomic(  51,  'Sb',  '[Cd]',  '5p-(2)5p(1)',           121.d0, 4.6802d0), &
       atomic(  52,  'Te',  '[Cd]',  '5p-(2)5p(2)',           130.d0, 4.7423d0), &
       atomic(  53,  'I ',  '[Cd]',  '5p-(2)5p(3)',           127.d0, 4.7500d0), &
       atomic(  54,  'Xe',  '[Cd]',  '5p-(2)5p(4)',           136.d0, 4.7964d0), &
       atomic(  55,  'Cs',  '[Xe]',  '6s(1)',                 133.d0, 4.8041d0), &
       atomic(  56,  'Ba',  '[Xe]',  '6s(2)',                 138.d0, 4.8378d0), &
       atomic(  57,  'La',  '[Xe]',  '5d-(1)6s(2)',           139.d0, 4.8550d0), &
       atomic(  58,  'Ce',  '[Xe]',  '4f-(1)5d-(1)6s(2)',     140.d0, 4.8771d0), &
       atomic(  59,  'Pr',  '[Xe]',  '4f-(3)6s(2)',           141.d0, 4.8919d0), &
       atomic(  60,  'Nd',  '[Xe]',  '4f-(4)6s(2)',           142.d0, 4.9123d0), &
       atomic(  61,  'Pm',  '[Xe]',  '4f-(5)6s(2)',           145.d0, 0.0000d0), &  !  ???
       atomic(  62,  'Sm',  '[Xe]',  '4f-(6)6s(2)',           144.d0, 4.9524d0), &
       atomic(  63,  'Eu',  '[Xe]',  '4f-(6)4f(1)6s(2)',      145.d0, 4.9663d0), &
       atomic(  64,  'Gd',  '[Xe]',  '4f-(6)4f(1)5d-(1)6s(2)',160.d0, 5.1734d0), &
       atomic(  65,  'Tb',  '[Xe]',  '4f-(6)4f(3)6s(2)',      159.d0, 5.0600d0), &
       atomic(  66,  'Dy',  '[Xe]',  '4f-(6)4f(4)6s(2)',      148.d0, 5.0455d0), &
       atomic(  67,  'Ho',  '[Xe]',  '4f-(6)4f(5)6s(2)',      165.d0, 5.2022d0), &
       atomic(  68,  'Er',  '[Xe]',  '4f-(6)4f(6)6s(2)',      170.d0, 5.2789d0), &
       atomic(  69,  'Tm',  '[Xe]',  '4f-(5)4f(8)6s(2)',      169.d0, 5.2256d0), &
       atomic(  70,  'Yb',  '[Xe]',  '4f-(6)4f(8)6s(2)',      176.d0, 5.3215d0), &
       atomic(  71,  'Lu',  '[4f]',  '5d-(1)6s(2)',           175.d0, 5.3700d0), &
       atomic(  72,  'Hf',  '[4f]',  '5d-(2)6s(2)',           178.d0, 5.3371d0), &
       atomic(  73,  'Ta',  '[4f]',  '5d-(3)6s(2)',           181.d0, 5.3507d0), &
       atomic(  74,  'W ',  '[4f]',  '5d-(4)6s(2)',           184.d0, 5.3658d0), &
       atomic(  75,  'Re',  '[4f]',  '5d-(4)5d(1)6s(2)',      185.d0, 5.3596d0), &
       atomic(  76,  'Os',  '[4f]',  '5d-(4)5d(2)6s(2)',      192.d0, 5.4126d0), &
       atomic(  77,  'Ir',  '[4f]',  '5d-(4)5d(3)6s(2)',      191.d0, 5.3968d0), &
       atomic(  78,  'Pt',  '[4f]',  '5d-(4)5d(5)6s(1)',      194.d0, 5.4236d0), &
       atomic(  79,  'Au',  '[4f]',  '5d-(4)5d(6)6s(1)',      197.d0, 5.4371d0), &
       atomic(  80,  'Hg',  '[4f]',  '5d-(4)5d(6)6s(2)',      198.d0, 5.4463d0), &
       atomic(  81,  'Tl',  '[Hg]',  '6p-(1)',                205.d0, 5.4759d0), &
       atomic(  82,  'Pb',  '[Hg]',  '6p-(2)',                208.d0, 5.5012d0), &
       atomic(  83,  'Bi',  '[Hg]',  '6p-(2)6p(1)',           209.d0, 5.5211d0), &
       atomic(  84,  'Po',  '[Hg]',  '6p-(2)6p(2)',           208.d0, 5.5584d0), &
       atomic(  85,  'At',  '[Hg]',  '6p-(2)6p(3)',           212.d0, 0.0000d0), &  ! ???
       atomic(  86,  'Rn',  '[Hg]',  '6p-(2)6p(4)',           212.d0, 5.5915d0), &
       atomic(  87,  'Fr',  '[Rn]',  '7s(1)',                 212.d0, 5.5915d0), &
       atomic(  88,  'Ra',  '[Rn]',  '7s(2)',                 214.d0, 5.6079d0), &
       atomic(  89,  'Ac',  '[Rn]',  '6d-(1)7s(2)',           227.d0, 0.0000d0), &  ! ???
       atomic(  90,  'Th',  '[Rn]',  '6d-(2)7s(2)',           232.d0, 5.7848d0), &
       atomic(  91,  'Pa',  '[Rn]',  '5f-(2)6d-(1)7s(2)',     231.d0, 0.0000d0), &  ! ???
       atomic(  92,  'U ',  '[Rn]',  '5f-(3)6d-(1)7s(2)',     238.d0, 5.8571d0), &
       atomic(  93,  'Np',  '[Rn]',  '5f-(4)6d-(1)7s(2)',     237.d0, 0.0000d0), &  ! ???
       atomic(  94,  'Pu',  '[Rn]',  '5f-(6)7s(2)',           239.d0, 5.8601d0), &
       atomic(  95,  'Am',  '[Rn]',  '5f-(6)5f(1)7s(2)',      243.d0, 5.9048d0), &
       atomic(  96,  'Cm',  '[Rn]',  '5f-(6)5f(1)6d-(1)7s(2)',244.d0, 5.8429d0), &
       atomic(  97,  'Bk',  '[Rn]',  '5f-(6)5f(3)7s(2)',      247.d0, 0.0000d0), &  ! ???  5.8160  where I took them?
       atomic(  98,  'Cf',  '[Rn]',  '5f-(6)5f(4)7s(2)',      251.d0, 0.0000d0), &  !      5.8440
       atomic(  99,  'Es',  '[Rn]',  '5f-(6)5f(5)7s(2)',      254.d0, 0.0000d0), &  !      5.8650
       atomic( 100,  'Fm',  '[Rn]',  '5f-(6)5f(6)7s(2)',      257.d0, 0.0000d0), &  !      5.8860
       atomic( 101,  'Md',  '[Rn]',  '5f-(6)5f(7)7s(2)',      258.d0, 0.0000d0), &  !      5.8930
       atomic( 102,  'No',  '[Rn]',  '5f-(6)5f(8)7s(2)',      259.d0, 0.0000d0), &  !      5.8860
       atomic( 103,  'Lr',  '[Rn]',  '5f-(6)5f(8)7s(2)7p-(1)',260.d0, 0.0000d0), &  !      5.9060
       atomic( 104,  'Rf',  '[Rn]',  '5f-(6)5f(8)6d(-2)7s(2)',265.d0, 0.0000d0) /)  !      5.8720

      Character(100) :: &
       He='1s', &
       Be='1s 2s', &
       Ne='1s 2s 2p- 2p', &
       Mg='1s 2s 2p- 2p 3s', &
       Ar='1s 2s 2p- 2p 3s 3p- 3p', &
       Zn='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s', &
       Kr='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p', &
       Cd='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 5s', &
       Xe='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 5s 5p- 5p', &
       f4='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 4f- 4f 5s 5p- 5p', &
       Hg='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 4f- 4f 5s 5p- 5p 5d- 5d 6s', &
       Rn='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 4f- 4f 5s 5p- 5p 5d- 5d 6s 6p- 6p', &
       f5='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 4f- 4f 5s 5p- 5p 5d- 5d 5f- 5f 6s 6p- 6p'

      End Module atoms_par



!======================================================================
      Subroutine Def_atom(an,atom,atw,rms,core,conf)
!======================================================================
! ... define the ground state of atom
!----------------------------------------------------------------------
      Use atoms_par

      Integer :: an, i
      Real(8) :: atw, rms
      Character(*) :: atom,core,conf

! ... if atomic number = 0, first try to define atom by symbol:

      if(an.le.0.or.an.gt.n_atoms) then
       an = -1
       Do i = 1,n_atoms
        if(atom.ne.atoms(i)%symbol) Cycle
        an = i; Exit
       End do
       if(an.eq.-1) Return
      else
       atom = atoms(an)%symbol
      end if

! ... atomic weight:

      atw = atoms(an)%weight

! ... atomic radius:

      rms = atoms(an)%rrms
      if(rms.eq.0.d0) then
       if(an.le.90)  then
         rms=(0.836d0*atw**(1.d0/3.d0)+0.570d0)
       else
         rms=(0.77d0*atw**(1.d0/3.d0)+0.980d0)
       end if
      end if

! ... core label:

      core = atoms(an)%core

      if(core(1:1).eq.'[') then
       Select case(core(2:3))
         case ('He'); core = trim(He)
         case ('Be'); core = trim(Be)
         case ('Ne'); core = trim(Ne)
         case ('Mg'); core = trim(Mg)
         case ('Ar'); core = trim(Ar)
         case ('Zn'); core = trim(Zn)
         case ('Kr'); core = trim(Kr)
         case ('Cd'); core = trim(Cd)
         case ('Xe'); core = trim(Xe)
         case ('4f'); core = trim(f4)
         case ('Hg'); core = trim(Hg)
         case ('Rn'); core = trim(Rn)
         CASE DEFAULT
         Stop 'Def_atom: problem with core'
       End Select
      end if

! ... configuration label:

      conf = atoms(an)%conf

      End Subroutine Def_atom
