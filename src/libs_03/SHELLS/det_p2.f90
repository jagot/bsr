!=======================================================================
  Subroutine Det_p2 (id,ML,MS,Idet)
!=======================================================================
! determinamt id and corresponding ML,MS for subshell p2 
!-----------------------------------------------------------------------
 
  Implicit none
 
  Integer, intent(in)  :: id
  Integer, intent(out) :: ML,MS,Idet(*)
 
  Integer, parameter :: iq_p2 =   2
  Integer, parameter :: kd_p2 =  15
 
  Integer :: Idet_p2 (iq_p2,kd_p2)
 
  Integer :: ML_p2 (kd_p2)
  Integer :: MS_p2 (kd_p2)
 
  if(id.le.0.or.id.gt.kd_p2) Stop "Det_p2: index id is out of range"
 
  ML = ML_p2 (id)
  MS = MS_p2 (id)
 
  Idet (1:iq_p2)= Idet_p2 (:,id)
 

  Data Idet_p2 ( 1,:)/ &
   1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5/

  Data Idet_p2 ( 2,:)/ &
   2, 3, 4, 5, 6, 3, 4, 5, 6, 4, 5, 6, 5, 6, 6/

  Data ML_p2 / &
  1, -1, -1,  3,  3, -1, -1,  3,  3, -3,  1,  1,  1,  1,  5/

  Data MS_p2 / &
  1, -1,  1, -1,  1,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1/
 
  End Subroutine Det_p2 
