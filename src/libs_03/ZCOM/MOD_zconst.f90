!=======================================================================
      Module zconst
!=======================================================================
!     contains basic processor-dependent constants as well as
!     a larger variety of physical constants
!-----------------------------------------------------------------------

! ... kind-constants

      integer, parameter :: int4 = selected_int_kind(9)
      integer, parameter :: int2 = selected_int_kind(4)
      integer, parameter :: int1 = selected_int_kind(2)

      integer, parameter :: sp = kind(1.0)
      integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))
      integer, parameter :: qp = selected_real_kind(2*precision(1.0_dp))

! ... simple constants

      real(dp), parameter :: zero  = 0.0_dp,    &
                             one   = 1.0_dp,    &
                             two   = 2.0_dp,    &
                             three = 3.0_dp,    &
                             four  = 4.0_dp,    &
                             five  = 5.0_dp,    &
                             six   = 6.0_dp,    &
                             seven = 7.0_dp,    &
                             eight = 8.0_dp,    &
                             nine  = 9.0_dp,    &
                             ten   = 10.0_dp,   &
                             half  = one/two,   &
                             third = one/three

! ... some useful physical constants
! ... (based on NIST compilation 2001)

      real(dp), parameter ::                            &
       bohr_radius_in_cm      =    0.5291772083e-8_dp,  &
       ao_cm                  =    0.5291772083e-8_dp,  &
       one_over_alpha         =    137.03599976e+0_dp,  &
       alpha                  =    7.297352533e-3_dp,   &
       c_au                   =    137.03599976e+0_dp,  &
       c_vacuum_in_m_per_s    =    299792458e+0_dp,     &
       electron_charge_in_C   =    1.602176462e-19_dp,  &
       electron_mass_in_g     =    9.10938188e-28_dp,   &
       electron_mass_in_amu   =    5.485799110e-4_dp,   &
       proton_mass_in_amu     =    1.00727646688e+0_dp, &
       proton_electron        =    1836.1526675e+0_dp,  &
       electron_proton        =    5.446170232e-4_dp,   &
       hbar_in_J_s            =    1.054571596e-27_dp,  &
       hbar_in_eV_s           =    6.58211889e-16_dp,   &
       au_eV_inf              =    27.2113834e+0_dp,    &
       au_cm_inf              =    219471.62e+0_dp,     &
       Ry_eV_inf              =    13.60569175e+0_dp,   &
       Ry_cm_inf              =    109737.31568549_dp,  &
       time_au_in_sec         =    4.134138e+16_dp,     &
       time_au                =    2.4189e-17_dp,       &
       pi_a0_2                =    0.87973e+0_dp,       &
       fermi_in_cm            =    1.0e-13_dp,          &
       pi                     =    3.141592653589793238462643e+0_dp

! ... Standard input, output, error and scratch unit numbers

      integer(int4) :: istdi=5, istdo=6 ,istde=0

! ... maximum length for record:

      integer(int4), parameter ::  mrecl = 100000

! ... packing bases:

      Integer(int4), parameter :: ibf = 2**4  ! overlap factors
      Integer(int4), parameter :: ibd = 2**15 ! overlap determinants

      Integer(int4), parameter :: jb = 10     ! max. number of shells
      Integer(int4), parameter :: jb2=jb**2, jb4=jb**4, jb8=jb**8

      Integer(int4) :: ib1 = 2,     ib2 = 2**2,  ib3 = 2**3,  ib4 = 2**4, &
                       ib5 = 2**5,  ib6 = 2**6,  ib7 = 2**7,  ib8 = 2**8, &
                       ib9 = 2**9,  ib10= 2**10, ib15= 2**15, ib20= 2**20

      End Module zconst

