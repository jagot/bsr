!=======================================================================
  Subroutine Det_p3 (id,ML,MS,Idet)
!=======================================================================
! determinamt id and corresponding ML,MS for subshell p3 
!-----------------------------------------------------------------------
 
  Implicit none
 
  Integer, intent(in)  :: id
  Integer, intent(out) :: ML,MS,Idet(*)
 
  Integer, parameter :: iq_p3 =   3
  Integer, parameter :: kd_p3 =  20
 
  Integer :: Idet_p3 (iq_p3,kd_p3)
 
  Integer :: ML_p3 (kd_p3)
  Integer :: MS_p3 (kd_p3)
 
  if(id.le.0.or.id.gt.kd_p3) Stop "Det_p3: index id is out of range"
 
  ML = ML_p3 (id)
  MS = MS_p3 (id)
 
  Idet (1:iq_p3)= Idet_p3 (:,id)
 

  Data Idet_p3 ( 1,:)/ &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4/

  Data Idet_p3 ( 2,:)/ &
   2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 3, 3, 3, 4, 4, 5, 4, 4, 5, 5/

  Data Idet_p3 ( 3,:)/ &
   3, 4, 5, 6, 4, 5, 6, 5, 6, 6, 4, 5, 6, 5, 6, 6, 5, 6, 6, 6/

  Data ML_p3 / &
 -1, -1,  3,  3, -3,  1,  1,  1,  1,  5, -3,  1,  1,  1,  1,  5, -1, -1,  3,  3/

  Data MS_p3 / &
  0,  2,  0,  2,  0, -2,  0,  0,  2,  0,  2,  0,  2,  2,  4,  2,  0,  2,  0,  2/
 
  End Subroutine Det_p3 
