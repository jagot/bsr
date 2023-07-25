!=======================================================================
  Real(8) Function gen_f3 (ip,id) 
!=======================================================================
!
! contains f3-shell cfp coefficients as  sqrt(m1/m2) 
!
! ip - index of parent term
! id - index of daughter term
!
!-----------------------------------------------------------------------
 
  Implicit none
 
  Integer, parameter :: nd =   17  ! number of daughter terms
  Integer, parameter :: np =    7  ! number of parent terms
  Integer, parameter :: ng =   64  ! number of non-zero cfp's
 
  Integer(2) :: ind(0:nd)
  Integer(1) :: ipa(ng)
  Integer(4) :: m1(ng)
  Integer(4) :: m2(ng)
 
  Integer :: ip,id,i
 
  if(id.le.0) Stop "gen_f3: id < 0"
  if(ip.le.0) Stop "gen_f3: ip < 0"
  if(id.gt.nd) Stop "gen_f3: id > nd"
  if(ip.gt.np) Stop "gen_f3: ip > np"
 
  gen_f3 = 0.d0
  Do i=ind(id-1)+1,ind(id)
   if(ip.ne.ipa(i)) Cycle
   gen_f3 = sqrt(abs(DBLE(m1(i)))/DBLE(m2(i)))
   if(m1(i).lt.0) gen_f3 = -gen_f3
   Exit
  End Do
 
  Data ind/  &
    0,    1,    4,    7,   10,   12,   15,   20,   24,   31,   36,   42, &
   47,   51,   55,   59,   62,   64  /
 
 
  Data ipa/  &
  4,  2,  4,  6,  2,  4,  6,  2,  4,  6,  4,  6,  4,  3,  5,  2,  4,  6, &
  3,  5,  2,  6,  3,  5,  2,  4,  6,  1,  3,  5,  7,  2,  6,  3,  5,  7, &
  2,  4,  6,  3,  5,  7,  2,  6,  3,  5,  7,  4,  3,  5,  7,  6,  3,  5, &
  7,  4,  6,  5,  7,  6,  5,  7,  6,  7  /
 
 
  Data m1/  &
        1,       -3,       -2,      -22,        1,       -2,       11, &
       11,       -2,      -65,       -2,        7,       -1,        5, &
      -11,        3,       -7,       22,       -8,       33,       11, &
      -27,       33,        8,       -1,       -1,      -11,        2, &
       -5,       -1,      -13,       11,       -3,      -55,       -4, &
       91,      -11,       -7,       65,       55,     -125,       13, &
      -65,      -33,       13,       52,        5,       -1,       10, &
       65,      -91,       -1,       13,       -8,       -5,       -7, &
       -1,        3,        4,       -1,        8,      -17,       -1, &
        1  /
 
 
  Data m2/  &
        1,        7,        9,       63,       14,        3,       42, &
       42,        9,      126,        9,        9,        2,       21, &
       42,       49,       18,      441,       49,       98,       49, &
       98,       98,       49,       14,        6,       42,        7, &
      126,       14,      126,       28,       28,      252,       77, &
      396,      294,       18,      882,      294,     1078,       66, &
      196,      196,     1764,      147,       36,        2,      189, &
      462,      297,        2,       54,       33,      297,       18, &
        9,       22,       11,        2,       33,       66,        2, &
        2  /
 
  End Function gen_f3
