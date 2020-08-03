!====================================================================
    MODULE spline_param
!====================================================================
!   contains basic spline parameters
!--------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: grid_type = 1   !  grid's type 
    INTEGER :: nug = 99 
    Character(40) :: AF_grid = 'knot.dat'

    INTEGER :: ks = 8  ! order of B-splines
    INTEGER :: ns = 0  ! number of splines
    INTEGER :: nv = 0  ! number of intervals ( = ns-ks+1 )
    INTEGER :: ml = 0  ! number of intervals from 0 to 1 (=1/h)
    INTEGER :: me = 0  ! number of intervals in the exponential region

    REAL(8) :: h = 0.125    ! initial step in the knot sequence for z*r
    REAL(8) :: hmax = 1.00  ! maximum step, t(ns+1) - t(ns) 
    REAL(8) :: rmax = 40.00 ! border radius, t(ns+1)

    END MODULE spline_param



