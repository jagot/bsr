!=======================================================================
  Real(8) Function detc_d8 (id,it)
!=======================================================================
! coefficient for determinamt id and term it for subshell d8 
!-----------------------------------------------------------------------
 
  Implicit none
 
  Integer, intent(in)  :: id,it
 
  Integer, parameter :: kd_d8 =  45
  Integer, parameter :: nt_d8 =   5
 
  Integer :: INT_d8 (kd_d8,nt_d8)
 
  Integer :: Norm_d8  = 70
 
  if(id.le.0.or.id.gt.kd_d8) Stop "detc_d8: index id is out of range"
  if(it.le.0.or.it.gt.nt_d8) Stop "detc_d8: index it is out of range"
 
  detc_d8 = dfloat(INT_d8(id,it))/dfloat(Norm_d8)
 
  detc_d8 = dsqrt(dabs(detc_d8))
 
  if(INT_d8(id,it).lt.0) detc_d8=-detc_d8
 

  Data INT_d8 (:,   1)/ &
        0,        0,      -14,       14,        0,        0,        0,        0,        0,        0,        0,        0, &
        0,        0,        0,        0,        0,        0,        0,        0,       14,        0,        0,        0, &
        0,      -14,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0, &
        0,        0,        0,        0,        0,        0,        0,        0,       14  /

  Data INT_d8 (:,   2)/ &
        0,       56,       28,       28,       56,        0,        0,        0,       28,       14,        0,        0, &
       14,       28,        0,      -28,      -14,        0,        0,      -14,       -7,      -14,      -28,        0, &
        0,       -7,      -14,        0,        0,        0,        0,        0,       42,       21,      -42,      -21, &
        0,        0,        0,        0,       21,       42,      -21,      -42,        0  /

  Data INT_d8 (:,   3)/ &
        0,        0,      -20,       20,        0,        0,        0,        0,        0,       30,        0,        0, &
      -30,        0,      -30,        0,       30,        0,        0,        0,       -5,      -30,        0,        0, &
        0,        5,        0,      -30,        0,      -20,        0,      -20,        0,       -5,        0,       -5, &
       20,        0,       20,        0,        5,        0,        5,        0,      -20  /

  Data INT_d8 (:,   4)/ &
        0,       14,        7,        7,       14,        0,      -70,      -35,       42,       21,      -35,      -70, &
       21,       42,        0,      -42,      -21,       70,       35,       56,       28,      -21,      -42,       35, &
       70,       28,       56,        0,       70,       35,      -70,      -35,      -28,      -14,       28,       14, &
       35,       70,      -35,      -70,      -14,      -28,       14,       28,        0  /

  Data INT_d8 (:,   5)/ &
       70,        0,       -1,        1,        0,       70,        0,       35,        0,        5,      -35,        0, &
       -5,        0,       40,        0,        5,        0,       35,        0,      -16,       -5,        0,      -35, &
        0,       16,        0,       40,        0,      -15,        0,      -15,        0,       30,        0,       30, &
       15,        0,       15,        0,      -30,        0,      -30,        0,       36  /
 
  End Function detc_d8
