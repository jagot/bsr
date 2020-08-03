!=====================================================================
      Module param_LS
!=====================================================================
!     main parameters in library ZCONFLS
!---------------------------------------------------------------------
      Implicit none

      Integer, parameter :: msh = 25     ! max. number of shells behind core
      Integer, parameter :: ibc = 2**15  ! packing parameter
      Integer, parameter :: ksmax=61*61  ! max. set number,  see ASET in ELF4
       

!----------------------------------------------------------------------
!     the matrix elements under consideration:
!----------------------------------------------------------------------
!
!     noper     - number of different operators
!     ioper(:)  - pointer to required operators
!     joper(:)  - pointer to required operator for given configurations
!
!     IT_oper(:,:) - pointer on the done calculation for specific
!                    operators and given terms
!     JT_oper(:,:) - pointer on the required operators for given
!                    subset of term between two configurations
!
!     Operator(1)   -   overlaps
!     Operator(2)   -   kinatic energy
!     Operator(3)   -   two-electron electrostatic
!     Operator(4)   -   spin-orbit
!     Operator(5)   -   spin-other-orbit
!     Operator(6)   -   spin-spin
!     Operator(7)   -   orbit-orbit

!-----------------------------------------------------------------------

      Integer, parameter :: noper=7      
      Integer ioper(noper)/1,1,1,0,0,0,0/, joper(noper)
      Integer koper(noper)     !  MPI copy
      Real(8) :: coper(noper)

      Integer, allocatable :: JT_oper(:,:)
      Real(8), allocatable :: CT_oper(:,:)
      Integer, allocatable :: JD_oper(:,:)  ! MPI copy

      End Module param_LS


