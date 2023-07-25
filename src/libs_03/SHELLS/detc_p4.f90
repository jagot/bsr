!=======================================================================
  Real(8) Function detc_p4 (id,it)
!=======================================================================
! coefficient for determinamt id and term it for subshell p4 
!-----------------------------------------------------------------------
 
  Implicit none
 
  Integer, intent(in)  :: id,it
 
  Integer, parameter :: kd_p4 =  15
  Integer, parameter :: nt_p4 =   3
 
  Integer :: INT_p4 (kd_p4,nt_p4)
 
  Integer :: Norm_p4  = 6
 
  if(id.le.0.or.id.gt.kd_p4) Stop "detc_p4: index id is out of range"
  if(it.le.0.or.it.gt.nt_p4) Stop "detc_p4: index it is out of range"
 
  detc_p4 = dfloat(INT_p4(id,it))/dfloat(Norm_p4)
 
  detc_p4 = dsqrt(dabs(detc_p4))
 
  if(INT_p4(id,it).lt.0) detc_p4=-detc_p4
 

  Data INT_p4 (:,   1)/ &
        0,        0,       -2,        2,        0,        0,        0,        0,        0,        0,        0,        0, &
        0,        0,       -2/

  Data INT_p4 (:,   2)/ &
        0,        6,        3,        3,        6,        0,       -6,       -3,        6,        3,       -3,       -6, &
        3,        6,        0/

  Data INT_p4 (:,   3)/ &
        6,        0,       -1,        1,        0,        6,        0,        3,        0,        3,       -3,        0, &
       -3,        0,        4/
 
  End Function detc_p4
