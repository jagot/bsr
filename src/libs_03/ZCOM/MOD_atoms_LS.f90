!======================================================================
      Module atoms_LS
!======================================================================
!     atomic parameters according to periodic table
!     in LS notation
!----------------------------------------------------------------------
      Implicit none
      Integer, parameter :: n_atoms = 104
      Type atomic
       INTEGER :: an
       CHARACTER(2) :: symbol 
       CHARACTER(4) :: core 
       CHARACTER(40) :: conf 
      End type atomic

      Type(atomic), dimension(n_atoms), parameter, public :: atoms = (/ &
       atomic(   1,  'H ',  '    ',  '1s(1)'                 ), & 
       atomic(   2,  'He',  '    ',  '1s(2)'                 ), & 
       atomic(   3,  'Li',  '[He]',  '2s(1)'                 ), & 
       atomic(   4,  'Be',  '[He]',  '2s(2)'                 ), & 
       atomic(   5,  'B ',  '[He]',  '2s(2)2p(1)'            ), & 
       atomic(   6,  'C ',  '[He]',  '2s(2)2p(2)'            ), & 
       atomic(   7,  'N ',  '[He]',  '2s(2)2p(3)'            ), & 
       atomic(   8,  'O ',  '[He]',  '2s(2)2p(4)'            ), & 
       atomic(   9,  'F ',  '[He]',  '2s(2)2p(5)'            ), & 
       atomic(  10,  'Ne',  '[He]',  '2s(2)2p(6)'            ), & 
       atomic(  11,  'Na',  '[Ne]',  '3s(1)'                 ), & 
       atomic(  12,  'Mg',  '[Ne]',  '3s(2)'                 ), & 
       atomic(  13,  'Al',  '[Ne]',  '3s(2)3p(1)'            ), & 
       atomic(  14,  'Si',  '[Ne]',  '3s(2)3p(2)'            ), & 
       atomic(  15,  'P ',  '[Ne]',  '3s(2)3p(3)'            ), & 
       atomic(  16,  'S ',  '[Ne]',  '3s(2)3p(4)'            ), & 
       atomic(  17,  'Cl',  '[Ne]',  '3s(2)3p(5)'            ), & 
       atomic(  18,  'Ar',  '[Ne]',  '3s(2)3p(6)'            ), & 
       atomic(  19,  'K ',  '[Ar]',  '4s(1)'                 ), & 
       atomic(  20,  'Ca',  '[Ar]',  '4s(2)'                 ), & 
       atomic(  21,  'Sc',  '[Ar]',  '3d(1)4s(2)'            ), & 
       atomic(  22,  'Ti',  '[Ar]',  '3d(2)4s(2)'            ), & 
       atomic(  23,  'V ',  '[Ar]',  '3d(3)4s(2)'            ), & 
       atomic(  24,  'Cr',  '[Ar]',  '3d(5)4s(1)'            ), & 
       atomic(  25,  'Mn',  '[Ar]',  '3d(5)4s(2)'            ), & 
       atomic(  26,  'Fe',  '[Ar]',  '3d(6)4s(2)'            ), & 
       atomic(  27,  'Co',  '[Ar]',  '3d(9)4s(2)'            ), & 
       atomic(  28,  'Ni',  '[Ar]',  '3d(8)4s(2)'            ), & 
       atomic(  29,  'Cu',  '[Ar]',  '3d(10)4s(1)'           ), & 
       atomic(  30,  'Zn',  '[Ar]',  '3d(10)4s(2)'           ), & 
       atomic(  31,  'Ga',  '[Ar]',  '3d(10)4s(2)4p(1)'      ), & 
       atomic(  32,  'Ge',  '[Ar]',  '3d(10)4s(2)4p(2)'      ), & 
       atomic(  33,  'As',  '[Ar]',  '3d(10)4s(2)4p(3)'      ), & 
       atomic(  34,  'Se',  '[Ar]',  '3d(10)4s(2)4p(4)'      ), & 
       atomic(  35,  'Br',  '[Ar]',  '3d(10)4s(2)4p(5)'      ), & 
       atomic(  36,  'Kr',  '[Ar]',  '3d(10)4s(2)4p(6)'      ), & 
       atomic(  37,  'Rb',  '[Kr]',  '5s(1)'                 ), & 
       atomic(  38,  'Sr',  '[Kr]',  '5s(2)'                 ), & 
       atomic(  39,  'Y ',  '[Kr]',  '4d(1)5s(2)'            ), & 
       atomic(  40,  'Zr',  '[Kr]',  '4d(2)5s(2)'            ), & 
       atomic(  41,  'Nb',  '[Kr]',  '4d(4)5s(1)'            ), & 
       atomic(  42,  'Mo',  '[Kr]',  '4d(5)5s(1)'            ), & 
       atomic(  43,  'Tc',  '[Kr]',  '4d(5)5s(2)'            ), &  
       atomic(  44,  'Ru',  '[Kr]',  '4d(7)5s(1)'            ), & 
       atomic(  45,  'Rh',  '[Kr]',  '4d(8)5s(1)'            ), & 
       atomic(  46,  'Pd',  '[Kr]',  '4d(10)'                ), & 
       atomic(  47,  'Ag',  '[Kr]',  '4d(10)5s(1)'           ), & 
       atomic(  48,  'Cd',  '[Kr]',  '4d(10)5s(2)'           ), & 
       atomic(  49,  'In',  '[Kr]',  '4d(10)5s(2)5p(1)'      ), & 
       atomic(  50,  'Sn',  '[Kr]',  '4d(10)5s(2)5p(2)'      ), & 
       atomic(  51,  'Sb',  '[Kr]',  '4d(10)5s(2)5p(3)'      ), & 
       atomic(  52,  'Te',  '[Kr]',  '4d(10)5s(2)5p(4)'      ), & 
       atomic(  53,  'I ',  '[Kr]',  '4d(10)5s(2)5p(5)'      ), & 
       atomic(  54,  'Xe',  '[Kr]',  '4d(10)5s(2)5p(6)'      ), & 
       atomic(  55,  'Cs',  '[Xe]',  '6s(1)'                 ), & 
       atomic(  56,  'Ba',  '[Xe]',  '6s(2)'                 ), & 
       atomic(  57,  'La',  '[Xe]',  '5d(1)6s(2)'            ), & 
       atomic(  58,  'Ce',  '[Xe]',  '4f(1)5d(1)6s(2)'       ), & 
       atomic(  59,  'Pr',  '[Xe]',  '4f(3)6s(2)'            ), & 
       atomic(  60,  'Nd',  '[Xe]',  '4f(4)6s(2)'            ), & 
       atomic(  61,  'Pm',  '[Xe]',  '4f(5)6s(2)'            ), &  
       atomic(  62,  'Sm',  '[Xe]',  '4f(6)6s(2)'            ), &  
       atomic(  63,  'Eu',  '[Xe]',  '4f(7)6s(2)'            ), &  
       atomic(  64,  'Gd',  '[Xe]',  '4f(7)5d(1)6s(2)'       ), &  
       atomic(  65,  'Tb',  '[Xe]',  '4f(9)6s(2)'            ), &  
       atomic(  66,  'Dy',  '[Xe]',  '4f(10)6s(2)'           ), &  
       atomic(  67,  'Ho',  '[Xe]',  '4f(11)6s(2)'           ), &  
       atomic(  68,  'Er',  '[Xe]',  '4f(12)6s(2)'           ), &  
       atomic(  69,  'Tm',  '[Xe]',  '4f(13)6s(2)'           ), &  
       atomic(  70,  'Yb',  '[Xe]',  '4f(14)6s(2)'           ), &  
       atomic(  71,  'Lu',  '[4f]',  '5d(1)6s(2)'            ), &  
       atomic(  72,  'Hf',  '[4f]',  '5d(2)6s(2)'            ), &  
       atomic(  73,  'Ta',  '[4f]',  '5d(3)6s(2)'            ), &  
       atomic(  74,  'W ',  '[4f]',  '5d(4)6s(2)'            ), &  
       atomic(  75,  'Re',  '[4f]',  '5d(5)6s(2)'            ), &  
       atomic(  76,  'Os',  '[4f]',  '5d(6)6s(2)'            ), &  
       atomic(  77,  'Ir',  '[4f]',  '5d(7)6s(2)'            ), &  
       atomic(  78,  'Pt',  '[4f]',  '5d(9)6s(1)'            ), &  
       atomic(  79,  'Au',  '[4f]',  '5d(10)6s(1)'           ), &  
       atomic(  80,  'Hg',  '[4f]',  '5d(10)6s(2)'           ), &  
       atomic(  81,  'Tl',  '[4f]',  '5d(10)6s(2)6p(1)'      ), &  
       atomic(  82,  'Pb',  '[4f]',  '5d(10)6s(2)6p(2)'      ), &  
       atomic(  83,  'Bi',  '[4f]',  '5d(10)6s(2)6p(3)'      ), &  
       atomic(  84,  'Po',  '[4f]',  '5d(10)6s(2)6p(4)'      ), &  
       atomic(  85,  'At',  '[4f]',  '5d(10)6s(2)6p(5)'      ), &              
       atomic(  86,  'Rn',  '[4f]',  '5d(10)6s(2)6p(6)'      ), & 
       atomic(  87,  'Fr',  '[Rn]',  '7s(1)'                 ), & 
       atomic(  88,  'Ra',  '[Rn]',  '7s(2)'                 ), & 
       atomic(  89,  'Ac',  '[Rn]',  '6d(1)7s(2)'            ), &              
       atomic(  90,  'Th',  '[Rn]',  '6d(2)7s(2)'            ), &              
       atomic(  91,  'Pa',  '[Rn]',  '5f(2)6d(1)7s(2)'       ), &              
       atomic(  92,  'U ',  '[Rn]',  '5f(3)6d(1)7s(2)'       ), &              
       atomic(  93,  'Np',  '[Rn]',  '5f(4)6d(1)7s(2)'       ), &              
       atomic(  94,  'Pu',  '[Rn]',  '5f(6)7s(2)'            ), & 
       atomic(  95,  'Am',  '[Rn]',  '5f(7)7s(2)'            ), & 
       atomic(  96,  'Cm',  '[Rn]',  '5f(7)6d(1)7s(2)'       ), &  
       atomic(  97,  'Bk',  '[Rn]',  '5f(9)7s(2)'            ), &  
       atomic(  98,  'Cf',  '[Rn]',  '5f(10)7s(2)'           ), &  
       atomic(  99,  'Es',  '[Rn]',  '5f(11)7s(2)'           ), &  
       atomic( 100,  'Fm',  '[Rn]',  '5f(12)7s(2)'           ), &  
       atomic( 101,  'Md',  '[Rn]',  '5f(13)7s(2)'           ), &  
       atomic( 102,  'No',  '[Rn]',  '5f(14)7s(2)'           ), &  
       atomic( 103,  'Lr',  '[Rn]',  '5f(14)7s(2)7p(1)'      ), &  
       atomic( 104,  'Rf',  '[Rn]',  '5f(14)6d(2)7s(2)'      ) /)  
                                                                                 
      Character(100) :: &                                                              
       He='1s', &                                                                      
       s1='1s', &                                                                      
       Be='1s 2s', &                                                                    
       s2='1s 2s', &                                                                    
       Ne='1s 2s 2p', &                                                             
       p2='1s 2s 2p', &                                                             
       Mg='1s 2s 2p 3s', &                                                          
       s3='1s 2s 2p 3s', &                                                          
       Ar='1s 2s 2p 3s 3p', &                                                   
       p3='1s 2s 2p 3s 3p', &                                                   
       d3='1s 2s 2p 3s 3p 3d', &                                                                                   
       Zn='1s 2s 2p 3s 3p 3d 4s', &                                                                                   
       s4='1s 2s 2p 3s 3p 3d 4s', &                                                                                   
       Kr='1s 2s 2p 3s 3p 3d 4s 4p', &                                  
       p4='1s 2s 2p 3s 3p 3d 4s 4p', &                                  
       d4='1s 2s 2p 3s 3p 3d 4s 4p 4d', &                        
       Cd='1s 2s 2p 3s 3p 3d 4s 4p 4d 5s', &                        
       s5='1s 2s 2p 3s 3p 3d 4s 4p 4d 5s', &                        
       p5='1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p', &                 
       Xe='1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p', &                 
       f4='1s 2s 2p 3s 3p 3d 4s 4p 4d 4f 5s 5p', &          
       Hg='1s 2s 2p 3s 3p 3d 4s 4p 4d 4f 5s 5p 5d 6s', & 
       s6='1s 2s 2p 3s 3p 3d 4s 4p 4d 4f 5s 5p 5d 6s', & 
       Rn='1s 2s 2p 3s 3p 3d 4s 4p 4d 4f 5s 5p 5d 6s 6p', &
       p6='1s 2s 2p 3s 3p 3d 4s 4p 4d 4f 5s 5p 5d 6s 6p', &
       f5='1s 2s 2p 3s 3p 3d 4s 4p 4d 4f 5s 5p 5d 5f 6s 6p'
                                                                                            
      End Module atoms_LS
   


!======================================================================
      Subroutine Def_atom_LS(an,atom,core,conf)
!======================================================================
! ... define the ground state of atom 
!----------------------------------------------------------------------
      Use atoms_LS

      Integer :: an, i
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
  
! ... core label: 
 
      core = atoms(an)%core

      if(core(1:1).eq.'[') then
       Select case(core(2:3))
         case ('He','s1'); core = trim(He)  
         case ('Be','s2'); core = trim(Be)  
         case ('Ne','p2'); core = trim(Ne)  
         case ('Mg','s3'); core = trim(Mg)  
         case ('Ar','p3'); core = trim(Ar)  
         case (     'd3'); core = trim(d3)  
         case ('Zn','s4'); core = trim(Zn)  
         case ('Kr','p4'); core = trim(Kr)  
         case (     'd4'); core = trim(d4)  
         case ('Cd','s5'); core = trim(Cd)  
         case ('Xe','p5'); core = trim(Xe)  
         case (     '4f'); core = trim(f4)  
         case ('Hg','s6'); core = trim(Hg)  
         case ('Rn','p6'); core = trim(Rn)  
         case (     'f5'); core = trim(f5)  
         case default; Stop 'Def_atom_LS: problem with core'
       End Select       
      end if

! ... configuration label: 

      conf = atoms(an)%conf

      End Subroutine Def_atom_LS

