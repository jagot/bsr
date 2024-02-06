!=======================================================================
  Subroutine Det_d7 (id,ML,MS,Idet)
!=======================================================================
! determinamt id and corresponding ML,MS for subshell d7
!-----------------------------------------------------------------------

  Implicit none

  Integer, intent(in)  :: id
  Integer, intent(out) :: ML,MS,Idet(*)

  Integer, parameter :: iq_d7 =   7
  Integer, parameter :: kd_d7 = 120

  Integer :: Idet_d7 (iq_d7,kd_d7)

  Integer :: ML_d7 (kd_d7)
  Integer :: MS_d7 (kd_d7)

  if(id.le.0.or.id.gt.kd_d7) Stop "Det_d7: index id is out of range"

  ML = ML_d7 (id)
  MS = MS_d7 (id)

  Idet (1:iq_d7)= Idet_d7 (:,id)


  Data Idet_d7 ( 1,:)/ &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4  /

  Data Idet_d7 ( 2,:)/ &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, &
   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 5, 5  /

  Data Idet_d7 ( 3,:)/ &
   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, &
   4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 5, 5, 5, 5, 5, 6, 6, &
   4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 5, 5, 5, 5, 5, 6, 6, 5, 5, 5, 5, 5, 6, 6, 6  /

  Data Idet_d7 ( 4,:)/ &
   4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 5, 5, 5, 5, 5, 5, 5, &
   5, 5, 5, 6, 6, 6, 6, 7, 6, 6, 6, 6, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 6, 6, 6, 6, 7, 7, 6, 6, 6, 6, 7, 7, 7, &
   5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 6, 6, 6, 6, 7, 7, 6, 6, 6, 6, 7, 7, 7, 6, 6, 6, 6, 7, 7, 7, 7  /

  Data Idet_d7 ( 5,:)/ &
   5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 7, 7, 7, 8, 8, 6, 6, 6, 6, 6, 6, 7, &
   7, 7, 8, 7, 7, 7, 8, 8, 7, 7, 7, 8, 8, 8, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 7, 7, 7, 8, 8, 7, 7, 7, 8, 8, 8, 7, 7, 7, 8, 8, 8, 8, &
   6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 7, 7, 7, 8, 8, 7, 7, 7, 8, 8, 8, 7, 7, 7, 8, 8, 8, 8, 7, 7, 7, 8, 8, 8, 8, 8  /

  Data Idet_d7 ( 6,:)/ &
   6, 6, 6, 6, 7, 7, 7, 8, 8, 9, 7, 7, 7, 8, 8, 9, 8, 8, 9, 9, 7, 7, 7, 8, 8, 9, 8, 8, 9, 9, 8, 8, 9, 9, 9, 7, 7, 7, 8, 8, 9, 8, &
   8, 9, 9, 8, 8, 9, 9, 9, 8, 8, 9, 9, 9, 9, 7, 7, 7, 8, 8, 9, 8, 8, 9, 9, 8, 8, 9, 9, 9, 8, 8, 9, 9, 9, 9, 8, 8, 9, 9, 9, 9, 9, &
   7, 7, 7, 8, 8, 9, 8, 8, 9, 9, 8, 8, 9, 9, 9, 8, 8, 9, 9, 9, 9, 8, 8, 9, 9, 9, 9, 9, 8, 8, 9, 9, 9, 9, 9, 9  /

  Data Idet_d7 ( 7,:)/ &
   7, 8, 9,10, 8, 9,10, 9,10,10, 8, 9,10, 9,10,10, 9,10,10,10, 8, 9,10, 9,10,10, 9,10,10,10, 9,10,10,10,10, 8, 9,10, 9,10,10, 9, &
  10,10,10, 9,10,10,10,10, 9,10,10,10,10,10, 8, 9,10, 9,10,10, 9,10,10,10, 9,10,10,10,10, 9,10,10,10,10,10, 9,10,10,10,10,10,10, &
   8, 9,10, 9,10,10, 9,10,10,10, 9,10,10,10,10, 9,10,10,10,10,10, 9,10,10,10,10,10,10, 9,10,10,10,10,10,10,10  /

  Data ML_d7 / &
 -3, -3,  5,  5, -9, -1, -1, -1, -1,  7, -9, -1, -1, -1, -1,  7, -7, -7,  1,  1, -5,  3,  3,  3,  3, 11, -3, -3,  5,  5, -3, -3, &
  5,  5, -1, -5,  3,  3,  3,  3, 11, -3, -3,  5,  5, -3, -3,  5,  5, -1,  1,  1,  9,  9,  3,  3, -7,  1,  1,  1,  1,  9, -5, -5, &
  3,  3, -5, -5,  3,  3, -3, -1, -1,  7,  7,  1,  1, -1, -1,  7,  7,  1,  1,  5, -7,  1,  1,  1,  1,  9, -5, -5,  3,  3, -5, -5, &
  3,  3, -3, -1, -1,  7,  7,  1,  1, -1, -1,  7,  7,  1,  1,  5, -3, -3,  5,  5, -1, -1,  3,  3  /

  Data MS_d7 / &
  0,  2,  0,  2,  0, -2,  0,  0,  2,  0,  2,  0,  2,  2,  4,  2,  0,  2,  0,  2,  0, -2,  0,  0,  2,  0, -2,  0, -2,  0,  0,  2, &
  0,  2,  0,  2,  0,  2,  2,  4,  2,  0,  2,  0,  2,  2,  4,  2,  4,  2,  0,  2,  0,  2,  0,  2,  0, -2,  0,  0,  2,  0, -2,  0, &
 -2,  0,  0,  2,  0,  2,  0, -2,  0, -2,  0, -2,  0,  0,  2,  0,  2,  0,  2,  0,  2,  0,  2,  2,  4,  2,  0,  2,  0,  2,  2,  4, &
  2,  4,  2,  0,  2,  0,  2,  0,  2,  2,  4,  2,  4,  2,  4,  2,  0,  2,  0,  2,  0,  2,  0,  2  /

  End Subroutine Det_d7
