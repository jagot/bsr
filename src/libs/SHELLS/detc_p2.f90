!=======================================================================
  Real(8) Function detc_p2 (id,it)
!=======================================================================
! coefficient for determinamt id and term it for subshell p2
!-----------------------------------------------------------------------

  Implicit none

  Integer, intent(in)  :: id,it

  Integer, parameter :: kd_p2 =  15
  Integer, parameter :: nt_p2 =   3

  Integer :: INT_p2 (kd_p2,nt_p2)

  Integer :: Norm_p2  = 6

  if(id.le.0.or.id.gt.kd_p2) Stop "detc_p2: index id is out of range"
  if(it.le.0.or.it.gt.nt_p2) Stop "detc_p2: index it is out of range"

  detc_p2 = dfloat(INT_p2(id,it))/dfloat(Norm_p2)

  detc_p2 = dsqrt(dabs(detc_p2))

  if(INT_p2(id,it).lt.0) detc_p2=-detc_p2


  Data INT_p2 (:,   1)/ &
        2,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,       -2, &
        2,        0,        0/

  Data INT_p2 (:,   2)/ &
        0,        6,        3,       -6,       -3,        3,        6,       -3,       -6,        0,       -6,       -3, &
       -3,       -6,        0/

  Data INT_p2 (:,   3)/ &
       -4,        0,       -3,        0,       -3,        3,        0,        3,        0,       -6,        0,       -1, &
        1,        0,       -6/

  End Function detc_p2
