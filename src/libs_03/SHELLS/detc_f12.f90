!=======================================================================
  Real(8) Function detc_f12 (id,it)
!=======================================================================
! coefficient for determinamt id and term it for subshell f12
!-----------------------------------------------------------------------

  Implicit none

  Integer, intent(in)  :: id,it

  Integer, parameter :: kd_f12 =  91
  Integer, parameter :: nt_f12 =   7

  Integer :: INT_f12(kd_f12,nt_f12)

  Integer :: Norm_f12 = 924

  if(id.le.0.or.id.gt.kd_f12) Stop "detc_f12: index id is out of range"
  if(it.le.0.or.it.gt.nt_f12) Stop "detc_f12: index it is out of range"

  detc_f12 = dfloat(INT_f12(id,it))/dfloat(Norm_f12)

  detc_f12 = dsqrt(dabs(detc_f12))

  if(INT_f12(id,it).lt.0) detc_f12=-detc_f12


  Data INT_f12(:,   1)/ &
        0,        0,     -132,      132,        0,        0,        0,        0,        0,        0,        0,        0, &
        0,        0,        0,        0,        0,        0,        0,        0,      132,        0,        0,        0, &
        0,     -132,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0, &
        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0, &
        0,        0,        0,        0,        0,        0,     -132,        0,        0,        0,        0,        0, &
        0,        0,        0,      132,        0,        0,        0,        0,        0,        0,        0,        0, &
        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0, &
        0,        0,        0,        0,        0,        0,     -132  /

  Data INT_f12(:,   2)/ &
        0,     -594,     -297,     -297,     -594,        0,        0,        0,     -198,      -99,        0,        0, &
      -99,     -198,        0,      198,       99,        0,        0,      264,      132,       99,      198,        0, &
        0,      132,      264,        0,        0,        0,        0,        0,        0,        0,      330,      165, &
        0,        0,        0,        0,        0,        0,      165,      330,        0,        0,        0,        0, &
        0,     -330,     -165,        0,        0,      -66,      -33,        0,        0,        0,        0,     -165, &
     -330,        0,        0,      -33,      -66,        0,        0,        0,        0,        0,        0,        0, &
        0,        0,      396,      198,     -396,     -198,        0,        0,        0,        0,        0,        0, &
        0,        0,      198,      396,     -198,     -396,        0  /

  Data INT_f12(:,   3)/ &
        0,        0,      275,     -275,        0,        0,        0,        0,        0,     -275,        0,        0, &
      275,        0,        0,        0,     -275,        0,        0,        0,        0,      275,        0,        0, &
        0,        0,        0,        0,        0,        0,        0,      110,        0,        0,        0,      165, &
        0,        0,     -110,        0,        0,        0,     -165,        0,     -264,        0,      110,        0, &
        0,        0,      165,        0,        0,        0,      -99,     -110,        0,        0,        0,     -165, &
        0,        0,        0,       99,        0,     -264,        0,        0,        0,        0,        0,     -220, &
        0,     -220,        0,      -22,        0,      -22,        0,        0,        0,        0,      220,        0, &
      220,        0,       22,        0,       22,        0,     -176  /

  Data INT_f12(:,   4)/ &
        0,     -308,     -154,     -154,     -308,        0,        0,        0,     -616,     -308,        0,        0, &
     -308,     -616,        0,      616,      308,        0,        0,     -308,     -154,      308,      616,        0, &
        0,     -154,     -308,        0,        0,        0,      616,      308,     -616,     -308,        0,        0, &
        0,        0,      308,      616,     -308,     -616,        0,        0,        0,     -616,     -308,        0, &
        0,        0,        0,      616,      308,      308,      154,     -308,     -616,        0,        0,        0, &
        0,      308,      616,      154,      308,        0,      308,      154,     -308,     -154,      308,      154, &
     -308,     -154,     -308,     -154,      308,      154,      154,      308,     -154,     -308,      154,      308, &
     -154,     -308,     -154,     -308,      154,      308,        0  /

  Data INT_f12(:,   5)/ &
        0,        0,       54,      -54,        0,        0,        0,        0,        0,     -180,        0,        0, &
      180,        0,      420,        0,     -180,        0,        0,        0,      294,      180,        0,        0, &
        0,     -294,        0,      420,        0,      252,        0,      324,        0,       84,        0,     -192, &
     -252,        0,     -324,        0,      -84,        0,      192,        0,      240,        0,      324,        0, &
      252,        0,     -192,        0,       84,        0,        6,     -324,        0,     -252,        0,      192, &
        0,      -84,        0,       -6,        0,      240,        0,     -378,        0,     -378,        0,       18, &
        0,       18,        0,       90,        0,       90,      378,        0,      378,        0,      -18,        0, &
      -18,        0,      -90,        0,      -90,        0,      216  /

  Data INT_f12(:,   6)/ &
        0,      -22,      -11,      -11,      -22,        0,      924,      462,     -110,      -55,      462,      924, &
      -55,     -110,        0,      110,       55,     -924,     -462,     -352,     -176,       55,      110,     -462, &
     -924,     -176,     -352,        0,     -924,     -462,      308,      154,      308,      154,     -594,     -297, &
     -462,     -924,      154,      308,      154,      308,     -297,     -594,        0,     -308,     -154,      924, &
      462,      594,      297,     -308,     -154,     -550,     -275,     -154,     -308,      462,      924,      297, &
      594,     -154,     -308,     -275,     -550,        0,      616,      308,     -616,     -308,     -616,     -308, &
      616,      308,      220,      110,     -220,     -110,      308,      616,     -308,     -616,     -308,     -616, &
      308,      616,      110,      220,     -110,     -220,        0  /

  Data INT_f12(:,   7)/ &
     -924,        0,        1,       -1,        0,     -924,        0,     -462,        0,       -7,      462,        0, &
        7,        0,     -504,        0,       -7,        0,     -462,        0,       36,        7,        0,      462, &
        0,      -36,        0,     -504,        0,      210,        0,       28,        0,     -378,        0,     -105, &
     -210,        0,      -28,        0,      378,        0,      105,        0,     -420,        0,       28,        0, &
      210,        0,     -105,        0,     -378,        0,      225,      -28,        0,     -210,        0,      105, &
        0,      378,        0,     -225,        0,     -420,        0,      -84,        0,      -84,        0,      224, &
        0,      224,        0,     -350,        0,     -350,       84,        0,       84,        0,     -224,        0, &
     -224,        0,      350,        0,      350,        0,     -400  /

  End Function detc_f12
