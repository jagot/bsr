!--------------------------------------------------------------------
      Real(8) Function ZGEN (Q,L,N1,N2)
!--------------------------------------------------------------------
!     Coefficients of fractional parentage:   (l,q,v'S'L'|}l,q+1,vLS)
!     N1,N2 - index of v'L'S' and vLS terms in the corresponding lists
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: Q,L,N1,N2
      Integer IV(375), IPV(7)
      Real(8) :: R,S
      Real(8), external :: gen_p3,gen_p4,gen_p5
      Real(8), external :: gen_d3,gen_d4,gen_d5,gen_d6,gen_d7,gen_d8,gen_d9
      Real(8), external :: gen_f3,gen_f4,gen_f5,gen_f6,gen_f7
      Integer :: i,n,IA1,IL1,IS1,IA2,IL2,IS2
      Integer, external :: Iterm_LS

      ZGEN=1.d0

      IF(Q.EQ.0.or.Q.EQ.1.or.Q.EQ.4*L+1) RETURN
      if(Q.LT.0.OR.Q.GT.4*L+1) then
       write(*,'(a,i2,a,i2)') 'Stop in ZGEN   L=',L,'  Q=',Q
       Stop ' ' 
      end if

      Select case(L)
      Case(1)
        Select case(q+1)
         Case(3);   ZGEN = gen_p3(n1,n2)
         Case(4);   ZGEN = gen_p4(n1,n2)
         Case(5);   ZGEN = gen_p5(n1,n2)
        End Select
      Case(2)   
        Select case(q+1)
         Case(3);   ZGEN = gen_d3(n1,n2)
         Case(4);   ZGEN = gen_d4(n1,n2)
         Case(5);   ZGEN = gen_d5(n1,n2)
         Case(6);   ZGEN = gen_d6(n1,n2)
         Case(7);   ZGEN = gen_d7(n1,n2)
         Case(8);   ZGEN = gen_d8(n1,n2)
         Case(9);   ZGEN = gen_d9(n1,n2)
        End Select
      Case(3)
        if(q.ge.2*L+1) then
         i = Iterm_LS (L,Q,N1,IA1,IL1,IS1)
         i=4*L+2-q; if(i.ge.3) then; i=IPV(i); IA1=IV(i+n1); end if
         i = Iterm_LS (l,Q+1,N2,IA2,IL2,IS2)
         i=4*L+2-q-1; if(i.ge.3) then; i=IPV(i); IA2=IV(i+n2); end if
         i = IL1+IL2+IS1+IS2-4-L-L-IA1-IA2 

         if(mod(i,2).ne.0) then
          write(*,'(10i6)') q,  N1,IA1,IL1,IS1 
          write(*,'(10i6)') q+1,N2,IA2,IL2,IS2  
          Stop 'ZGEN: problems with phase'
         end if

         S =  (4*L+2-Q)*IL1*IS1
         S =  S/((Q+1)*IL2*IS2)
         S = sqrt(S) * (-1)**(i/2)
        end if
        Select case(q+1)
         Case(3);   ZGEN = gen_f3(n1,n2)
         Case(4);   ZGEN = gen_f4(n1,n2)
         Case(5);   ZGEN = gen_f5(n1,n2)
         Case(6);   ZGEN = gen_f6(n1,n2)
         Case(7);   ZGEN = gen_f7(n1,n2)
         Case(8);   ZGEN = gen_f7(n2,n1)*S 
         Case(9);   ZGEN = gen_f6(n2,n1)*S
         Case(10);  ZGEN = gen_f5(n2,n1)*S
         Case(11);  ZGEN = gen_f4(n2,n1)*S
         Case(12);  ZGEN = gen_f3(n2,n1)*S
         Case(13);  ZGEN = S
        End Select
      Case default
        write(*,'(a,i2,a,i2)') 'Stop in ZGEN:  L=',L,'  Q=',Q
        Stop ' ' 
      End Select

DATA IV /                                     &
 3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, &
 3, 3,                                        &
 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, &
 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 2, &
 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, &
 4, 4,                                        &
 5, 5, 5, 3, 5, 5, 3, 5, 5, 3, 5, 5, 5, 3, 5, &
 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 3, 5, 5, &
 5, 3, 3, 5, 5, 5, 1, 3, 5, 5, 5, 5, 5, 3, 3, &
 5, 5, 5, 5, 3, 3, 5, 5, 5, 5, 5, 3, 5, 5, 5, &
 5, 3, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5,       &
 6, 4, 6, 4, 6, 6, 4, 6, 4, 6, 6, 6, 6, 4, 6, &
 6, 6, 2, 4, 4, 6, 6, 6, 4, 4, 6, 6, 6, 2, 4, &
 4, 4, 6, 6, 6, 6, 6, 4, 4, 4, 6, 6, 6, 6, 2, &
 4, 4, 4, 6, 6, 6, 6, 6, 4, 4, 6, 6, 6, 6, 4, &
 4, 6, 6, 6, 6, 4, 6, 6, 4, 6, 6, 6, 6, 0, 4, &
 6, 6, 6, 2, 4, 4, 4, 6, 6, 4, 6, 6, 6, 2, 4, &
 4, 4, 6, 6, 6, 6, 4, 4, 6, 6, 2, 4, 4, 6, 6, &
 6, 6, 4, 6, 6, 4, 4, 6, 6, 6, 6, 4, 6, 6,    &
 7, 5, 7, 5, 7, 5, 7, 3, 7, 5, 5, 3, 5, 5, 7, &
 7, 7, 3, 5, 5, 5, 7, 3, 5, 5, 5, 7, 7, 7, 5, &
 5, 5, 7, 7, 3, 5, 5, 7, 7, 5, 5, 7, 5, 7, 7, &
 5, 7, 7, 7, 3, 5, 5, 5, 7, 3, 3, 5, 5, 5, 7, &
 7, 1, 3, 5, 5, 5, 5, 5, 7, 7, 7, 3, 3, 5, 5, &
 5, 5, 7, 7, 7, 7, 3, 3, 5, 5, 5, 5, 5, 7, 7, &
 3, 5, 5, 5, 5, 7, 7, 7, 7, 3, 5, 5, 5, 5, 7, &
 7, 3, 5, 5, 7, 7, 5, 5, 7, 7, 5, 7, 5, 7     /

 DATA IPV / 0,0,0, 17, 64, 137, 256 / 

 End Function ZGEN

!=======================================================================
  Real(8) Function gen_p3 (ip,id) 
!=======================================================================
!
! contains p3-shell cfp coefficients as  sqrt(m1/m2) 
!
! ip - index of parent term
! id - index of daughter term
!
!-----------------------------------------------------------------------

  Implicit none
  Integer, intent(in) :: ip,id
  Integer, parameter :: nd =   3  ! number of daughter terms
  Integer, parameter :: np =   3  ! number of parent terms
  Integer :: IN(nd),IC(np,nd)

  if(id.le.0) Stop "gen_p3: id < 0"
  if(ip.le.0) Stop "gen_p3: ip < 0"
  if(id.gt.nd) Stop "gen_p3: id > nd"
  if(ip.gt.np) Stop "gen_p3: ip > np"

  gen_p3=SQRT(abs(DBLE(IC(ip,id)))/DBLE(IN(id)))
  if(IC(ip,id).lt.0) gen_p3=-gen_p3

  Data IN/1,18,2/, IC/0,1,0,  4,-9,-5,  0,1,-1/

  End Function gen_p3
  

!=======================================================================
  Real(8) Function gen_p4 (ip,id) 
!=======================================================================
!
! contains p4-shell cfp coefficients as  sqrt(m1/m2) 
!
! ip - index of parent term
! id - index of daughter term
!
!-----------------------------------------------------------------------

  Implicit none
  Integer, intent(in) :: ip,id
  Integer, parameter :: nd =   3  ! number of daughter terms
  Integer, parameter :: np =   3  ! number of parent terms
  Integer :: IN(nd),IC(np,nd)

  if(id.le.0)  Stop "gen_p4: id < 0"
  if(ip.le.0)  Stop "gen_p4: ip < 0"
  if(id.gt.nd) Stop "gen_p4: id > nd"
  if(ip.gt.np) Stop "gen_p4: ip > np"

  gen_p4=SQRT(abs(DBLE(IC(ip,id)))/DBLE(IN(id)))
  if(IC(ip,id).lt.0) gen_p4=-gen_p4

  Data IN/1,12,4/, IC/0,1,0, -4,-3, 5,  0,-1,-3/

  End Function gen_p4


!=======================================================================
  Real(8) Function gen_p5 (ip,id) 
!=======================================================================
!
! contains p5-shell cfp coefficients as  sqrt(m1/m2) 
!
! ip - index of parent term
! id - index of daughter term
!
!-----------------------------------------------------------------------

  Implicit none
  Integer, intent(in) :: ip,id
  Integer, parameter :: nd =  1  ! number of daughter terms
  Integer, parameter :: np =  3  ! number of parent terms
  Integer :: IN(nd),IC(np,nd)

  if(id.le.0)  Stop "gen_p5: id < 0"
  if(ip.le.0)  Stop "gen_p5: ip < 0"
  if(id.gt.nd) Stop "gen_p5: id > nd"
  if(ip.gt.np) Stop "gen_p5: ip > np"

  gen_p5=SQRT(abs(DBLE(IC(ip,id)))/DBLE(IN(id)))
  if(IC(ip,id).lt.0) gen_p5=-gen_p5

  Data IN/15/, IC/1,9,5/

  End Function gen_p5


!=======================================================================
  Real(8) Function gen_d3 (ip,id) 
!=======================================================================
!
! contains d3-shell cfp coefficients as  sqrt(m1/m2) 
!
! ip - index of parent term
! id - index of daughter term
!
!-----------------------------------------------------------------------

  Implicit none
  Integer, intent(in) :: ip,id
  Integer, parameter :: nd =   8  ! number of daughter terms
  Integer, parameter :: np =   5  ! number of parent terms
  Integer :: IN(nd),IC(np,nd)

  if(id.le.0)  Stop "gen_d3: id < 0"
  if(ip.le.0)  Stop "gen_d3: ip < 0"
  if(id.gt.nd) Stop "gen_d3: id > nd"
  if(ip.gt.np) Stop "gen_d3: ip > np"

  gen_d3=SQRT(abs(DBLE(IC(ip,id)))/DBLE(IN(id)))
  if(IC(ip,id).lt.0) gen_d3=-gen_d3

  Data IN/30,15,60,140,70,5,42,2/

  Data IC/   0,  7, 15, -8,  0,   0, -8,  0, -7,  0, &
            16, -9, -5,-21, -9,   0,-49, 45, 21,-25, &
             0, 28,-10,  7,-25,   0, -1,  0,  4,  0, &
             0,  0,-10, 21, 11,   0,  0,  0, -1,  1  /

  End Function gen_d3


!=======================================================================
  Real(8) Function gen_d4 (ip,id) 
!=======================================================================
!
! contains d4-shell cfp coefficients as  sqrt(i1/i2) 
!
! ip - index of parent term
! id - index of daughter term
!
!-----------------------------------------------------------------------

  Implicit none
  Integer, intent(in) :: ip,id
  Integer, parameter :: nd =  16  ! number of daughter terms
  Integer, parameter :: np =   8  ! number of parent terms
  Integer :: IN(nd),IC(np,nd)

  if(id.le.0)  Stop "gen_d4: id < 0"
  if(ip.le.0)  Stop "gen_d4: ip < 0"
  if(id.gt.nd) Stop "gen_d4: id > nd"
  if(ip.gt.np) Stop "gen_d4: ip > np"

  gen_d4=SQRT(abs(DBLE(IC(ip,id)))/DBLE(IN(id)))
  if(IC(ip,id).lt.0) gen_d4=-gen_d4

  Data IN/    1,    1,  360,   90,   280,  140,  210,  10, &
            560,  840, 1680,  504,  1008, 1680,   60,  10/

  Data IC/    0,   0,   1,   0,   0,   0,   0,   0, &
              0,   0,   0,   1,   0,   0,   0,   0, &
            -14, -64, 135, -35, -56, -56,   0,   0, &
             25, -14,   0,  10, -25,  16,   0,   0, &
            -42,   0, 105,  45,  28,   0, -60,   0, &
             42,   0,   0,  20,  63,   0,  15,   0, &
            -14,  49,   0,  60, -21, -21,  45,   0, &
              0,   3,   0,   0,   0,   7,   0,   0, &
            120,   0,   0, 200,-105,   0,  -3,-132, &
             16, -56, 315,  15, -14, 224,  90, 110, &
           -200,-448,   0, 120,-175,-112,-405, 220, &
              0,   0, 189,- 25,  70,   0,  66,-154, &
              0,   0,   0,  88, 385,   0,-507, -28, &
              0,   0,   0, 200, 315,-560, 297, 308, &
              0,   0,   0,   0,   5,  20,  -9,  26, &
              0,   0,   0,   0,   0,   0,   3,   7  /

  End Function gen_d4



!=======================================================================
  Real(8) Function gen_d5 (ip,id) 
!=======================================================================
!
! contains d5-shell cfp coefficients as  sqrt(i1/i2) 
!
! ip - index of parent term
! id - index of daughter term
!
!-----------------------------------------------------------------------

  Implicit none
  Integer, intent(in) :: ip,id
  Integer, parameter :: nd =  16  ! number of daughter terms
  Integer, parameter :: np =  16  ! number of parent terms
  Integer :: IN(nd),IC(np,nd)

  if(id.le.0)  Stop "gen_d5: id < 0"
  if(ip.le.0)  Stop "gen_d5: ip < 0"
  if(id.gt.nd) Stop "gen_d5: id > nd"
  if(ip.gt.np) Stop "gen_d5: ip > np"

  gen_d5=SQRT(abs(DBLE(IC(ip,id)))/DBLE(IN(id)))
  if(IC(ip,id).lt.0) gen_d5=-gen_d5

  Data IN/   5,    1,  150,  300,    50,  350,  700, 700, &
          2800, 2800,  700, 8400, 18480,  420, 1100, 550/

  Data IC /5*0,-2,3,9*0,7*0,1,8*0,0,0,14,25,30,15,10,0,-15,      &
           -16,-25,5*0,0,0,-64,14,0,0,35,-75,0,-56,56,5*0,       &
           6,0,-9,0,-5,4*0,-21,0,-9,4*0,0,-14,-49,-14,45,-10,    &
           60,0,35,21,-21,-25,-11,45,0,0,0,-56,0,126,0,90,60,    &
           0,35,0,189,0,99,45,5*0,126,0,0,-135,-175,0,0,-84,     &
           0,0,180,4*0,448,-200,-160,180,120,0,105,112,-175,     &
          -400,275,-405,220,4*0,360,0,-100,600,0,-525,0,-315,    &
           0,495,-9,-396,0,0,0,-56,-16,0,0,-15,-175,0,224,14,    &
           0,0,-90,-110,0,4*0,-800,-100,600,0,-7,1680,945,880,   &
           845,891,924,-728,5*0,1452,968,0,-2541,0,4235,0,-1215, &
           -5577,-308,-2184,6*0,25,-105,0,0,-70,0,0,-66,154,0,   &
           8*0,33,-220,55,220,-5,-99,286,182,12*0,-45,99,231,-175/


  End Function gen_d5


!=======================================================================
  Real(8) Function gen_d6 (ip,id) 
!=======================================================================
!
! contains d6-shell cfp coefficients as  sqrt(i1/i2) 
!
! ip - index of parent term
! id - index of daughter term
!
!-----------------------------------------------------------------------

  Implicit none
  Integer, intent(in) :: ip,id
  Integer, parameter :: nd =  16  ! number of daughter terms
  Integer, parameter :: np =  16  ! number of parent terms
  Integer :: IN(nd),IC(np,nd)

  if(id.le.0)  Stop "gen_d6: id < 0"
  if(ip.le.0)  Stop "gen_d6: ip < 0"
  if(id.gt.nd) Stop "gen_d6: id > nd"
  if(ip.gt.np) Stop "gen_d6: ip > np"

  gen_d6=SQRT(abs(DBLE(IC(ip,id)))/DBLE(IN(id)))
  if(IC(ip,id).lt.0) gen_d6=-gen_d6

  Data IN/   1,    3,  270,  270,   210, 2100,  630,  30, &
          1680, 2520, 5040,  378, 33264, 5040, 1980, 330/

  Data IC /4*0,1,16*0,1,-2,9*0,0,0,14,64,-45,35,0,0,56,0,56,7*0, &
           25,-14,0,10,45,-90,-25,-45,16,7*0,42,0,-35,-45,0,0,   &
           -28,0,0,60,4*0,-280,0,210,0,0,100,450,0,315,175,0,75, &
           495,0,0,0,-42,0,-14,49,0,60,-30,-135,-21,105,-21,45,  &
           -33,75,0,0,0,6,0,3,0,0,0,5,0,0,7,0,0,9,4*0,120,0,0,   &
           200,-100,0,-105,-525,0,-3,495,0,-132,0,0,0,-64,224,   &
           -420,-60,0,0,56,0,-896,-360,0,0,-440,0,0,0,-200,-448, &
           0,120,540,480,-175,315,-112,-405,825,1200,220,0,4*0,  &
           -63,25,0,0,-70,0,0,-66,0,0,154,6*0,968,4356,0,4235,   &
           -7623,0,-5577,-3645,0,-308,-6552,5*0,200,-100,800,315,&
           -7,-560,297,845,-880,308,-728,8*0,55,99,220,-99,-15,  &
           -660,286,546,11*0,33,-45,0,77,-175/


  End Function gen_d6


!=======================================================================
  Real(8) Function gen_d7 (ip,id) 
!=======================================================================
!
! contains d7-shell cfp coefficients as  sqrt(i1/i2) 
!
! ip - index of parent term
! id - index of daughter term
!
!-----------------------------------------------------------------------

  Implicit none
  Integer, intent(in) :: ip,id
  Integer, parameter :: nd =   8  ! number of daughter terms
  Integer, parameter :: np =  16  ! number of parent terms
  Integer :: IN(nd),IC(np,nd)

  if(id.le.0)  Stop "gen_d7: id < 0"
  if(ip.le.0)  Stop "gen_d7: ip < 0"
  if(id.gt.nd) Stop "gen_d7: id > nd"
  if(ip.gt.np) Stop "gen_d7: ip > np"

  gen_d7=SQRT(abs(DBLE(IC(ip,id)))/DBLE(IN(id)))
  if(IC(ip,id).lt.0) gen_d7=-gen_d7

  Data IN /  210,  210,  140,  980,  1960,  490, 5880, 770/

  Data IC /  0,   0,   7, -50,  15, -30, -20,   0, &
            30,  -8,  50,   0,   0,   0,   0,   0, &
             0,   0, -16, -14,   0,   0, -35,  75, &
             0, -14, -56,   0,   0,   0,   0,   0, &
             8,   0,  27,   0,  15,   0,   0,   0, &
             0,  63,   0,  27,   0,   0,   0,   0, &
             0,  56, -49,  56,  45,  40,-240,   0, &
          -140,  21,  84, -25,  44,-180,   0,   0, &
             0,   0, 112, 200, -40,-180,-120,   0, &
          -105,  28, 175,-100,-275, 405,-220,   0, &
             0,   0, -14,  16,   0,   0,  15, 175, &
             0,  56, -14,   0,   0,  90, 110,   0, &
             0,   0,   0,   0,-200, 100,-600,   0, &
             7, 420,-945, 220,-845,-891,-924, 728, &
             0,   0,   0,   0,   0,   0,   0,   0, &
           -33, -55, -55,  55,   5,  99,-286,-182  /

  End Function gen_d7

!=======================================================================
  Real(8) Function gen_d8 (ip,id) 
!=======================================================================
!
! contains d8-shell cfp coefficients as  sqrt(i1/i2) 
!
! ip - index of parent term
! id - index of daughter term
!
!-----------------------------------------------------------------------

  Implicit none
  Integer, intent(in) :: ip,id
  Integer, parameter :: nd =   5  ! number of daughter terms
  Integer, parameter :: np =   8  ! number of parent terms
  Integer :: IN(nd),IC(np,nd)

  if(id.le.0)  Stop "gen_d8: id < 0"
  if(ip.le.0)  Stop "gen_d8: ip < 0"
  if(id.gt.nd) Stop "gen_d8: id > nd"
  if(ip.gt.np) Stop "gen_d8: ip > np"

  gen_d8=SQRT(abs(DBLE(IC(ip,id)))/DBLE(IN(id)))
  if(IC(ip,id).lt.0) gen_d8=-gen_d8

  Data IN /  1,  240,  560,  560, 336/

  Data IC /  0,   0,   1,   0,   0,   0,   0,   0, &
           -14, -64, -15, -35, -56, -56,   0,   0, &
          -126,   0, -35, 135,  84,   0,-180,   0, &
            16, -56, -35,  15, -14, 224,  90, 110, &
             0,   0, -21, -25,  70,   0,  66,-154  /

  End Function gen_d8


!=======================================================================
  Real(8) Function gen_d9 (ip,id) 
!=======================================================================
!
! contains d9-shell cfp coefficients as  sqrt(i1/i2) 
!
! ip - index of parent term
! id - index of daughter term
!
!-----------------------------------------------------------------------

  Implicit none
  Integer, intent(in) :: ip,id
  Integer, parameter :: nd =   1 ! number of daughter terms
  Integer, parameter :: np =   5 ! number of parent terms
  Integer :: IN(nd),IC(np,nd)

  if(id.le.0)  Stop "gen_d9: id < 0"
  if(ip.le.0)  Stop "gen_d9: ip < 0"
  if(id.gt.nd) Stop "gen_d9: id > nd"
  if(ip.gt.np) Stop "gen_d9: ip > np"

  gen_d9=SQRT(abs(DBLE(IC(ip,id)))/DBLE(IN(id)))
  if(IC(ip,id).lt.0) gen_d9=-gen_d9

  Data IN / 45 /

  Data IC / 1, 9, 5, 21, 9/

  End Function gen_d9
