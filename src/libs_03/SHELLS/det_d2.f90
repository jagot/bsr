!=======================================================================
  Subroutine Det_d2 (id,ML,MS,Idet)
!=======================================================================
! determinamt id and corresponding ML,MS for subshell d2 
!-----------------------------------------------------------------------
 
  Implicit none
 
  Integer, intent(in)  :: id
  Integer, intent(out) :: ML,MS,Idet(*)
 
  Integer, parameter :: iq_d2 =   2
  Integer, parameter :: kd_d2 =  45
 
  Integer :: Idet_d2 (iq_d2,kd_d2)
 
  Integer :: ML_d2 (kd_d2)
  Integer :: MS_d2 (kd_d2)
 
  if(id.le.0.or.id.gt.kd_d2) Stop "Det_d2: index id is out of range"
 
  ML = ML_d2 (id)
  MS = MS_d2 (id)
 
  Idet (1:iq_d2)= Idet_d2 (:,id)
 

  Data Idet_d2 ( 1,:)/ &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, &
   8, 8, 9  /

  Data Idet_d2 ( 2,:)/ &
   2, 3, 4, 5, 6, 7, 8, 9,10, 3, 4, 5, 6, 7, 8, 9,10, 4, 5, 6, 7, 8, 9,10, 5, 6, 7, 8, 9,10, 6, 7, 8, 9,10, 7, 8, 9,10, 8, 9,10, &
   9,10,10  /

  Data ML_d2 / &
  1, -1, -1,  3,  3, -3, -3,  5,  5, -1, -1,  3,  3, -3, -3,  5,  5, -3,  1,  1, -5, -5,  3,  3,  1,  1, -5, -5,  3,  3,  5, -1, &
 -1,  7,  7, -1, -1,  7,  7, -7,  1,  1,  1,  1,  9  /

  Data MS_d2 / &
  1, -1,  1, -1,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1, -1, &
  1, -1,  1,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1  /
 
  End Subroutine Det_d2 
