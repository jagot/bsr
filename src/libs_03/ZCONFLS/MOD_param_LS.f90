!=====================================================================
      Module param_LS
!=====================================================================
!     main parameters in library ZCONFLS
!---------------------------------------------------------------------
      Implicit none

      Integer, parameter :: msh = 25     ! max. number of shells behind core
      Integer, parameter :: ibc = 2**15  ! packing parameter
      Integer, parameter :: noper=7      ! different operators
      Integer, parameter :: ksmax=61*61  ! max. set number,  see ASET in ELF4

      End Module param_LS


