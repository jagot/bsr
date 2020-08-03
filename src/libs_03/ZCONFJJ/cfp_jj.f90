!=====================================================================
      Real(8) Function cfp_jj(j,q,JP,VP,JD,VD)
!=====================================================================
!     fractional parentage coefficients for j^q subshells
!     j  - 2*j value 
!     q  - the number of electrons in daughter subshell
!     JP - 2*J value of parent subshell
!     VP - seniority of parent term
!     JD - 2*J value of daughter subshell
!     VD - seniority of daughter term
!     coefficents are stored in form: +- sqrt(i1/i2)
!     subshell terms are defined in routine Jterm
!---------------------------------------------------------------------
      Implicit none

      Integer :: j,q,JP,VP,JD,VD, WP,QP, WD,QD, IP,ID
      Real(8) :: CN,CD,phase
      Real(8) :: zero = 0.d0, one = 1.d0
      Integer, external :: Jterm

      Integer, parameter :: n2_j3=2, n3_j3=1            
      Integer :: num3_j3(n3_j3,n2_j3), norm3_j3(n3_j3)

      Integer, parameter :: n2_j5=3, n3_j5=3            
      Integer :: num3_j5(n3_j5,n2_j5), norm3_j5(n3_j5)

      Integer, parameter :: n2_j7=4, n3_j7=6, n4_j7=8
      Integer :: num3_j7(n3_j7,n2_j7), norm3_j7(n3_j7)
      Integer :: num4_j7(n4_j7,n3_j7), norm4_j7(n4_j7)

! ... check of j and q values

      if(q.lt.1.or.q.gt.j+1) then
       write(*,*) 'cfp_jj: q is out of scope:',q 
       Stop 'Stop in cfp_jj'
      end if

      if(j.lt.1.or.(j.ge.9.and.q.gt.2)) then
       write(*,*) 'cfp_jj: j is out of scope:',j,q 
       Stop 'Stop in cfp_jj'
      end if

! ... term indexes:

      IP = Jterm (j,q-1,0,JP,VP,WP,QP)
      ID = Jterm (j,q,  0,JD,VD,WD,QD)

! .... trivial cases:

      if(q.le.2.or.q.eq.j+1) then
       cfp_jj = one
       Return
      end if

! ... select different subshells:

      phase = one 
      Select case(j*100+q)  
      Case(303)                                  !  3/2 ^ 3
       CN = num3_j3(ID,IP)
       CD = norm3_j3(ID)
       if(CN.lt.zero) phase=-phase
      Case(503)                                  !  5/2 ^ 3
       CN = num3_j5(ID,IP)
       CD = norm3_j5(ID)
       if(CN.lt.zero) phase=-phase
      Case(504)                                  !  5/2 ^ 4
       CN = num3_j5(IP,ID)
       CD = norm3_j5(IP)
       if(CN.lt.zero) phase=-phase
       CN = CN * ((7-q)*(1+JP))
       CD = CD * (q*(1+JD))
       phase = phase * (-1)**((JD-JP-VD+VP)/2-3)
      Case(505)                                  !  5/2 ^ 5
       CN = ((7-q)*(1+JP))
       CD = (q*(1+JD))
       phase =  (-1)**((JD-JP-VD+VP)/2-3)
      Case(703)                                  !  7/2 ^ 3
       CN = num3_j7(ID,IP)
       CD = norm3_j7(ID)
       if(CN.lt.zero) phase=-phase
      Case(704)                                  !  7/2 ^ 4
       CN = num4_j7(ID,IP)
       CD = norm4_j7(ID)
       if(CN.lt.zero) phase=-phase
      Case(705)                                  !  7/2 ^ 5
       CN = num4_j7(IP,ID)
       CD = norm4_j7(IP)
       if(CN.lt.zero) phase=-phase
       CN = CN * ((9-q)*(1+JP))
       CD = CD * (q*(1+JD))
       phase = phase * (-1)**((JD-JP-VD+VP)/2-3)
      Case(706)                                  !  7/2 ^ 6
       CN = num3_j7(IP,ID)
       CD = norm3_j7(IP)
       if(CN.lt.zero) phase=-phase
       CN = CN * ((9-q)*(1+JP))
       CD = CD * (q*(1+JD))
       phase = phase * (-1)**((JD-JP-VD+VP)/2-3)
      Case(707)                                  !  7/2 ^ 7
       CN = ((9-q)*(1+JP))
       CD = (q*(1+JD))
       phase =  (-1)**((JD-JP-VD+VP)/2-3)
      Case default
       write(*,*) 'cfp_jj: unknown subshell j^q: ',j,q
       Stop 'Stop in cfp_jj'
      End Select

      cfp_jj = sqrt(abs(CN)/CD) * phase
      Return

      DATA num3_j3 /  1, -5 /
      DATA norm3_j3/  6/

      DATA num3_j5 / -4, 0,  0, &
                      5,-5,  3, &
                      9, 2,-11 /
      DATA norm3_j5/ 18, 7, 14 /

      DATA num3_j7 /  9,   0,   0,   0,   0,   0, & 
                     -5,   3, 121, 143, -55,   0, & 
                     -9, -11,  12,-900,  39,   5, & 
                    -13,   0, -65, 343, 104, -17  / 
      DATA norm3_j7/ 36,  14, 198,1386, 198,  22  /

      DATA num4_j7 /  1,  280,  308, 1144,    0,    0,    0,   0, &
                      0,   54, -121,    0, -968,  169,  462,   0, &
                      0, -231,  -14,  195,  -77, 2366, -343,   0, &
                      0,  -65,  250, -245,-1755,   90, -945, 140, &
                      0, -210,   91,  624,  280, 2275,  650, 234, &
                      0,    0,  140,-1224,    0, -560,  680, 627  /
      DATA norm4_j7/  1,  840,  924, 3432, 3080, 5460, 3080,1001  /

      End Function cfp_jj

