!=======================================================================
  Subroutine Det_f4 (id,ML,MS,Idet)
!=======================================================================
! determinamt id and corresponding ML,MS for subshell f4
!-----------------------------------------------------------------------

  Implicit none

  Integer, intent(in)  :: id
  Integer, intent(out) :: ML,MS,Idet(*)

  Integer, parameter :: iq_f4 =   4
  Integer, parameter :: kd_f4 =1001

  Integer :: Idet_f4 (iq_f4,kd_f4)

  Integer :: ML_f4 (kd_f4)
  Integer :: MS_f4 (kd_f4)

  if(id.le.0.or.id.gt.kd_f4) Stop "Det_f4: index id is out of range"

  ML = ML_f4 (id)
  MS = MS_f4 (id)

  Idet (1:iq_f4)= Idet_f4 (:,id)


  Data Idet_f4 ( 1,:)/ &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
   2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, &
   4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
   4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
   4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, &
   5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
   5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, &
   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
   6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
   8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,11/

  Data Idet_f4 ( 2,:)/ &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, &
   4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, &
   5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, &
   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, &
   8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,11,11,11,12, 3, 3, 3, 3, 3, 3, 3, 3, &
   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
   3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
   4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
   5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
   7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,11,11, &
  11,12, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
   4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, &
   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
   7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,11,11,11,12, 5, &
   5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, &
   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
   8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,11,11,11,12, 6, 6, 6, 6, 6, 6, 6, &
   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
   8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,11,11,11,12, 7, 7, 7, 7, 7, 7, 7, &
   7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10, &
  10,10,10,11,11,11,12, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,11,11,11,12, &
   9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,11,11,11,12,10,10,10,10,10,10,11,11,11,12,11,11,11,12,12/

  Data Idet_f4 ( 3,:)/ &
   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, &
   7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, &
   5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 5, 5, 5, 5, 5, &
   5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 6, 6, &
   6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 7, 7, 7, 7, 7, 7, 7, 8, &
   8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 9, &
   9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13,10,10,10,10,11,11,11,12,12,13,11,11,11,12,12,13,12,12,13,13, 4, 4, 4, 4, 4, 4, 4, 4, &
   4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11, &
  11,11,12,12,13, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10, &
  10,10,11,11,11,12,12,13, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12, &
  12,13, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10, &
  10,10,10,11,11,11,12,12,13, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13,10,10,10,10,11,11,11,12,12,13,11,11,11,12,12,13,12,12, &
  13,13, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11, &
  11,11,12,12,13, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 7, &
   7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10, &
  11,11,11,12,12,13, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13,10,10,10,10,11,11,11,12,12,13,11,11,11,12,12,13,12,12,13,13, 6, &
   6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 7, 7, 7, 7, 7, 7, 7, &
   8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, &
   9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13,10,10,10,10,11,11,11,12,12,13,11,11,11,12,12,13,12,12,13,13, 7, 7, 7, 7, 7, 7, 7, &
   8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, &
   9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13,10,10,10,10,11,11,11,12,12,13,11,11,11,12,12,13,12,12,13,13, 8, 8, 8, 8, 8, 8, 9, &
   9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13,10,10,10,10,11,11,11,12,12,13,11,11,11, &
  12,12,13,12,12,13,13, 9, 9, 9, 9, 9,10,10,10,10,11,11,11,12,12,13,10,10,10,10,11,11,11,12,12,13,11,11,11,12,12,13,12,12,13,13, &
  10,10,10,10,11,11,11,12,12,13,11,11,11,12,12,13,12,12,13,13,11,11,11,12,12,13,12,12,13,13,12,12,13,13,13/

  Data Idet_f4 ( 4,:)/ &
   4, 5, 6, 7, 8, 9,10,11,12,13,14, 5, 6, 7, 8, 9,10,11,12,13,14, 6, 7, 8, 9,10,11,12,13,14, 7, 8, 9,10,11,12,13,14, 8, 9,10,11, &
  12,13,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, 5, 6, 7, 8, 9,10,11,12,13,14, 6, 7, 8, 9,10,11,12,13, &
  14, 7, 8, 9,10,11,12,13,14, 8, 9,10,11,12,13,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, 6, 7, 8, 9,10, &
  11,12,13,14, 7, 8, 9,10,11,12,13,14, 8, 9,10,11,12,13,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, 7, 8, &
   9,10,11,12,13,14, 8, 9,10,11,12,13,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, 8, 9,10,11,12,13,14, 9, &
  10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14,10, &
  11,12,13,14,11,12,13,14,12,13,14,13,14,14,11,12,13,14,12,13,14,13,14,14,12,13,14,13,14,14,13,14,14,14, 5, 6, 7, 8, 9,10,11,12, &
  13,14, 6, 7, 8, 9,10,11,12,13,14, 7, 8, 9,10,11,12,13,14, 8, 9,10,11,12,13,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12, &
  13,14,13,14,14, 6, 7, 8, 9,10,11,12,13,14, 7, 8, 9,10,11,12,13,14, 8, 9,10,11,12,13,14, 9,10,11,12,13,14,10,11,12,13,14,11,12, &
  13,14,12,13,14,13,14,14, 7, 8, 9,10,11,12,13,14, 8, 9,10,11,12,13,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13, &
  14,14, 8, 9,10,11,12,13,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, 9,10,11,12,13,14,10,11,12,13,14,11, &
  12,13,14,12,13,14,13,14,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14,11,12,13,14,12,13,14,13,14,14,12,13,14,13,14,14,13,14, &
  14,14, 6, 7, 8, 9,10,11,12,13,14, 7, 8, 9,10,11,12,13,14, 8, 9,10,11,12,13,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12, &
  13,14,13,14,14, 7, 8, 9,10,11,12,13,14, 8, 9,10,11,12,13,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, 8, &
   9,10,11,12,13,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14, &
  12,13,14,13,14,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14,11,12,13,14,12,13,14,13,14,14,12,13,14,13,14,14,13,14,14,14, 7, &
   8, 9,10,11,12,13,14, 8, 9,10,11,12,13,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, 8, 9,10,11,12,13,14, &
   9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, &
  10,11,12,13,14,11,12,13,14,12,13,14,13,14,14,11,12,13,14,12,13,14,13,14,14,12,13,14,13,14,14,13,14,14,14, 8, 9,10,11,12,13,14, &
   9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, 9,10,11,12,13,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14, &
  10,11,12,13,14,11,12,13,14,12,13,14,13,14,14,11,12,13,14,12,13,14,13,14,14,12,13,14,13,14,14,13,14,14,14, 9,10,11,12,13,14,10, &
  11,12,13,14,11,12,13,14,12,13,14,13,14,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14,11,12,13,14,12,13,14,13,14,14,12,13,14, &
  13,14,14,13,14,14,14,10,11,12,13,14,11,12,13,14,12,13,14,13,14,14,11,12,13,14,12,13,14,13,14,14,12,13,14,13,14,14,13,14,14,14, &
  11,12,13,14,12,13,14,13,14,14,12,13,14,13,14,14,13,14,14,14,12,13,14,13,14,14,13,14,14,14,13,14,14,14,14/

  Data ML_f4 / &
 -3,  1,  1, -5, -5,  3,  3, -7, -7,  5,  5,  1,  1, -5, -5,  3,  3, -7, -7,  5,  5,  5, -1, -1,  7,  7, -3, -3,  9,  9, -1, -1, &
  7,  7, -3, -3,  9,  9, -7,  1,  1, -9, -9,  3,  3,  1,  1, -9, -9,  3,  3,  9, -1, -1, 11, 11, -1, -1, 11, 11,-11,  1,  1,  1, &
  1, 13, -1, -1, -7, -7,  1,  1, -9, -9,  3,  3,  3, -3, -3,  5,  5, -5, -5,  7,  7, -3, -3,  5,  5, -5, -5,  7,  7, -9, -1, -1, &
-11,-11,  1,  1, -1, -1,-11,-11,  1,  1,  7, -3, -3,  9,  9, -3, -3,  9,  9,-13, -1, -1, -1, -1, 11,  3, -3, -3,  5,  5, -5, -5, &
  7,  7, -3, -3,  5,  5, -5, -5,  7,  7, -9, -1, -1,-11,-11,  1,  1, -1, -1,-11,-11,  1,  1,  7, -3, -3,  9,  9, -3, -3,  9,  9, &
-13, -1, -1, -1, -1, 11,  1,  1,  9,  9, -1, -1, 11, 11, -5,  3,  3, -7, -7,  5,  5,  3,  3, -7, -7,  5,  5, 11,  1,  1, 13, 13, &
  1,  1, 13, 13, -9,  3,  3,  3,  3, 15, -5,  3,  3, -7, -7,  5,  5,  3,  3, -7, -7,  5,  5, 11,  1,  1, 13, 13,  1,  1, 13, 13, &
 -9,  3,  3,  3,  3, 15, -3, -3,-13,-13, -1, -1,  5, -5, -5,  7,  7, -5, -5,  7,  7,-15, -3, -3, -3, -3,  9,  5, -5, -5,  7,  7, &
 -5, -5,  7,  7,-15, -3, -3, -3, -3,  9,  3,  3, 15, 15, -7,  5,  5,  5,  5, 17, -7,  5,  5,  5,  5, 17, -5, -5,  7,  7, -1, -1, &
 -7, -7,  1,  1, -9, -9,  3,  3,  3, -3, -3,  5,  5, -5, -5,  7,  7, -3, -3,  5,  5, -5, -5,  7,  7, -9, -1, -1,-11,-11,  1,  1, &
 -1, -1,-11,-11,  1,  1,  7, -3, -3,  9,  9, -3, -3,  9,  9,-13, -1, -1, -1, -1, 11,  3, -3, -3,  5,  5, -5, -5,  7,  7, -3, -3, &
  5,  5, -5, -5,  7,  7, -9, -1, -1,-11,-11,  1,  1, -1, -1,-11,-11,  1,  1,  7, -3, -3,  9,  9, -3, -3,  9,  9,-13, -1, -1, -1, &
 -1, 11,  1,  1,  9,  9, -1, -1, 11, 11, -5,  3,  3, -7, -7,  5,  5,  3,  3, -7, -7,  5,  5, 11,  1,  1, 13, 13,  1,  1, 13, 13, &
 -9,  3,  3,  3,  3, 15, -5,  3,  3, -7, -7,  5,  5,  3,  3, -7, -7,  5,  5, 11,  1,  1, 13, 13,  1,  1, 13, 13, -9,  3,  3,  3, &
  3, 15, -3, -3,-13,-13, -1, -1,  5, -5, -5,  7,  7, -5, -5,  7,  7,-15, -3, -3, -3, -3,  9,  5, -5, -5,  7,  7, -5, -5,  7,  7, &
-15, -3, -3, -3, -3,  9,  3,  3, 15, 15, -7,  5,  5,  5,  5, 17, -7,  5,  5,  5,  5, 17, -5, -5,  7,  7,  1, -5, -5,  3,  3, -7, &
 -7,  5,  5, -5, -5,  3,  3, -7, -7,  5,  5,-11, -3, -3,-13,-13, -1, -1, -3, -3,-13,-13, -1, -1,  5, -5, -5,  7,  7, -5, -5,  7, &
  7,-15, -3, -3, -3, -3,  9, -1, -1,  7,  7, -3, -3,  9,  9, -7,  1,  1, -9, -9,  3,  3,  1,  1, -9, -9,  3,  3,  9, -1, -1, 11, &
 11, -1, -1, 11, 11,-11,  1,  1,  1,  1, 13, -7,  1,  1, -9, -9,  3,  3,  1,  1, -9, -9,  3,  3,  9, -1, -1, 11, 11, -1, -1, 11, &
 11,-11,  1,  1,  1,  1, 13, -5, -5,-15,-15, -3, -3,  3, -7, -7,  5,  5, -7, -7,  5,  5,-17, -5, -5, -5, -5,  7,  3, -7, -7,  5, &
  5, -7, -7,  5,  5,-17, -5, -5, -5, -5,  7,  1,  1, 13, 13, -9,  3,  3,  3,  3, 15, -9,  3,  3,  3,  3, 15, -7, -7,  5,  5, -1, &
 -1,  7,  7, -3, -3,  9,  9, -7,  1,  1, -9, -9,  3,  3,  1,  1, -9, -9,  3,  3,  9, -1, -1, 11, 11, -1, -1, 11, 11,-11,  1,  1, &
  1,  1, 13, -7,  1,  1, -9, -9,  3,  3,  1,  1, -9, -9,  3,  3,  9, -1, -1, 11, 11, -1, -1, 11, 11,-11,  1,  1,  1,  1, 13, -5, &
 -5,-15,-15, -3, -3,  3, -7, -7,  5,  5, -7, -7,  5,  5,-17, -5, -5, -5, -5,  7,  3, -7, -7,  5,  5, -7, -7,  5,  5,-17, -5, -5, &
 -5, -5,  7,  1,  1, 13, 13, -9,  3,  3,  3,  3, 15, -9,  3,  3,  3,  3, 15, -7, -7,  5,  5, -3,  5,  5, -5, -5,  7,  7,  5,  5, &
 -5, -5,  7,  7, 13,  3,  3, 15, 15,  3,  3, 15, 15, -7,  5,  5,  5,  5, 17, -1, -1,-11,-11,  1,  1,  7, -3, -3,  9,  9, -3, -3, &
  9,  9,-13, -1, -1, -1, -1, 11,  7, -3, -3,  9,  9, -3, -3,  9,  9,-13, -1, -1, -1, -1, 11,  5,  5, 17, 17, -5,  7,  7,  7,  7, &
 19, -5,  7,  7,  7,  7, 19, -3, -3,  9,  9, -1, -1,-11,-11,  1,  1,  7, -3, -3,  9,  9, -3, -3,  9,  9,-13, -1, -1, -1, -1, 11, &
  7, -3, -3,  9,  9, -3, -3,  9,  9,-13, -1, -1, -1, -1, 11,  5,  5, 17, 17, -5,  7,  7,  7,  7, 19, -5,  7,  7,  7,  7, 19, -3, &
 -3,  9,  9,  1, -9, -9,  3,  3, -9, -9,  3,  3,-19, -7, -7, -7, -7,  5, -1, -1, 11, 11,-11,  1,  1,  1,  1, 13,-11,  1,  1,  1, &
  1, 13, -9, -9,  3,  3, -1, -1, 11, 11,-11,  1,  1,  1,  1, 13,-11,  1,  1,  1,  1, 13, -9, -9,  3,  3, -3,  9,  9,  9,  9, 21, &
 -1, -1, 11, 11, -1, -1, 11, 11,  1/

  Data MS_f4 / &
  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1, -1,  1, -1,  1,  1,  3, &
  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1, -1,  1,  1, &
  3,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1, -3, -1, -3, -1, -3, -1, -3, -1, -1,  1, -1,  1, -1,  1, -1,  1, -1, -3, -1, &
 -3, -1, -3, -1, -1,  1, -1,  1, -1,  1, -1, -3, -1, -3, -1, -1,  1, -1,  1, -1, -3, -1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, &
 -1,  1,  1,  3,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3, &
  1, -1,  1,  1,  3,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1, -3, -1, -3, -1, -3, -1, -1,  1, -1,  1, -1,  1, -1, -3, -1, -3, -1, &
 -1,  1, -1,  1, -1, -3, -1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3, &
  1, -1,  1,  1,  3,  1, -1,  1, -1,  1, -1,  1, -1, -3, -1, -3, -1, -1,  1, -1,  1, -1, -3, -1, -1,  1, -1,  1, -1,  1, -1,  1, &
  1,  3,  1,  3,  1, -1,  1,  1,  3,  1, -1,  1, -1,  1, -1, -3, -1, -1,  1, -1,  1, -1,  1,  1,  3,  1, -1,  1, -1,  1,  1,  3, &
  1,  3,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1, -1,  1, &
  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1,  3,  1,  3,  1,  3,  3,  5, &
  3,  5,  3,  5,  3,  5,  3,  1,  3,  1,  3,  1,  3,  3,  5,  3,  5,  3,  5,  3,  1,  3,  1,  3,  3,  5,  3,  5,  3,  1,  3,  3, &
  5,  3,  1,  3,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3, &
  1, -1,  1,  1,  3,  1,  3,  1,  3,  1,  3,  1,  3,  3,  5,  3,  5,  3,  5,  3,  1,  3,  1,  3,  3,  5,  3,  5,  3,  1,  3,  3, &
  5,  3,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1,  3,  3,  5,  3,  5, &
  3,  1,  3,  3,  5,  3,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1,  3,  1,  3,  3,  5,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1, -1, &
  1, -1,  1,  1,  3,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1, &
  3,  1, -1,  1,  1,  3,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1, -3, -1, -3, -1, -3, -1, -1,  1, -1,  1, -1,  1, -1, -3, -1, -3, &
 -1, -1,  1, -1,  1, -1, -3, -1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1, &
  3,  1, -1,  1,  1,  3,  1, -1,  1, -1,  1, -1,  1, -1, -3, -1, -3, -1, -1,  1, -1,  1, -1, -3, -1, -1,  1, -1,  1, -1,  1, -1, &
  1,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1, -1,  1, -1,  1, -1, -3, -1, -1,  1, -1,  1, -1,  1,  1,  3,  1, -1,  1, -1,  1,  1, &
  3,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1, -1,  1, &
  1,  3,  1,  3,  1,  3,  1,  3,  1,  3,  3,  5,  3,  5,  3,  5,  3,  1,  3,  1,  3,  3,  5,  3,  5,  3,  1,  3,  3,  5,  3,  1, &
  3,  1,  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1,  3,  3,  5,  3,  5,  3,  1,  3, &
  3,  5,  3,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1,  3,  1,  3,  3,  5,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1, -1,  1,  1,  3, &
  1,  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1, -1,  1, -1,  1, -1,  1, -1, -3, -1, -3, -1, -1,  1, &
 -1,  1, -1, -3, -1, -1,  1, -1,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1, -1,  1, -1,  1, -1, -3, -1, -1,  1, &
 -1,  1, -1,  1,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1, &
  3,  1,  3,  1,  3,  3,  5,  3,  5,  3,  1,  3,  3,  5,  3,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1,  3,  1,  3,  3,  5,  3,  1, &
  3,  1,  3,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1, -1,  1, -1,  1, -1, -3, -1, -1,  1, -1,  1, -1,  1,  1, &
  3,  1, -1,  1, -1,  1,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1,  3,  1,  3,  3,  5,  3,  1,  3,  1,  3,  1, -1,  1,  1,  3,  1, &
 -1,  1, -1,  1,  1,  3,  1,  3,  1/

  End Subroutine Det_f4
