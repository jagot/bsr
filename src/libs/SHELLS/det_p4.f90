!=======================================================================
  Subroutine Det_p4 (id,ML,MS,Idet)
!=======================================================================
! determinamt id and corresponding ML,MS for subshell p4
!-----------------------------------------------------------------------

  Implicit none

  Integer, intent(in)  :: id
  Integer, intent(out) :: ML,MS,Idet(*)

  Integer, parameter :: iq_p4 =   4
  Integer, parameter :: kd_p4 =  15

  Integer :: Idet_p4 (iq_p4,kd_p4)

  Integer :: ML_p4 (kd_p4)
  Integer :: MS_p4 (kd_p4)

  if(id.le.0.or.id.gt.kd_p4) Stop "Det_p4: index id is out of range"

  ML = ML_p4 (id)
  MS = MS_p4 (id)

  Idet (1:iq_p4)= Idet_p4 (:,id)


  Data Idet_p4 ( 1,:)/ &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3/

  Data Idet_p4 ( 2,:)/ &
   2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 3, 3, 3, 4, 4/

  Data Idet_p4 ( 3,:)/ &
   3, 3, 3, 4, 4, 5, 4, 4, 5, 5, 4, 4, 5, 5, 5/

  Data Idet_p4 ( 4,:)/ &
   4, 5, 6, 5, 6, 6, 5, 6, 6, 6, 5, 6, 6, 6, 6/

  Data ML_p4 / &
 -3,  1,  1,  1,  1,  5, -1, -1,  3,  3, -1, -1,  3,  3,  1/

  Data MS_p4 / &
  1, -1,  1,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1/

  End Subroutine Det_p4
