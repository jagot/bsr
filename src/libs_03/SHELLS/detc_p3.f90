!=======================================================================
  Real(8) Function detc_p3 (id,it)
!=======================================================================
! coefficient for determinamt id and term it for subshell p3
!-----------------------------------------------------------------------

  Implicit none

  Integer, intent(in)  :: id,it

  Integer, parameter :: kd_p3 =  20
  Integer, parameter :: nt_p3 =   3

  Integer :: INT_p3 (kd_p3,nt_p3)

  Integer :: Norm_p3  = 6

  if(id.le.0.or.id.gt.kd_p3) Stop "detc_p3: index id is out of range"
  if(it.le.0.or.it.gt.nt_p3) Stop "detc_p3: index it is out of range"

  detc_p3 = dfloat(INT_p3(id,it))/dfloat(Norm_p3)

  detc_p3 = dsqrt(dabs(detc_p3))

  if(INT_p3(id,it).lt.0) detc_p3=-detc_p3


  Data INT_p3 (:,   1)/ &
        0,        0,        0,        0,        0,        6,        2,        2,        2,        0,        0,        2, &
        2,        2,        6,        0,        0,        0,        0,        0/

  Data INT_p3 (:,   2)/ &
        3,        3,        3,        3,        0,        0,       -3,        3,        0,        0,        0,        0, &
       -3,        3,        0,        0,        3,        3,        3,        3/

  Data INT_p3 (:,   3)/ &
        3,        3,       -3,       -3,       -6,        0,       -1,       -1,       -4,        6,       -6,        4, &
        1,        1,        0,        6,       -3,       -3,        3,        3/

  End Function detc_p3
