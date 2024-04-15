!--------------------------------------------------------------------
!     Routines for spectroscopic labelling of one-electron orbitals:
!
!     AL  - spectroscopic symbol for given L
!     LA  - value L from spectroscopic symbol
!     ELF4 - specroscopic notation for electron orbital (n,l,k) (a4)
!     ELF3 - specroscopic notation for electron orbital (n,l,k) (a3)
!     EL4_NLK(EL,n,l,k) -  decodes the specroscopic notation for a4
!     EL3_NLK(EL,n,l,k) -  decodes the specroscopic notation for a3
!--------------------------------------------------------------------

!====================================================================
      Character Function AL(L,k)
!====================================================================
!     provides spectroscopic symbols for L values ( L <= 153 )
!     k=1 - small letters; k=2 - capital letters
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: L,K
      Integer :: I
      Character(21) :: AS, AB

      DATA AS/'spdfghiklmnoqrtuvwxyz'/
      DATA AB/'SPDFGHIKLMNOQRTUVWXYZ'/

      i = L+1; IF(k.eq.5.or.k.eq.6) i=(L-1)/2+1
      if(i.ge.1.and.i.le.21) then
       if(k.eq.1.or.k.eq.5) AL=AS(I:I)
       if(k.eq.2.or.k.eq.6) AL=AB(I:I)
      elseif(i.ge.22.and.i.le.153) then
       AL=CHAR(i+101)  ! from '{' and futher
      else
       write(*,*) 'L,k=',L,k
       Stop ' AL: L is out of range'
      end if

      End Function AL


!====================================================================
      Integer Function LA(a)
!====================================================================
!     gives the value of L from spetroscopic symbol "a"
!--------------------------------------------------------------------
      Implicit none
      Character, Intent(in) :: a
      Character(21) :: AS, AB
      Integer :: i

      DATA AS/'spdfghiklmnoqrtuvwxyz'/
      DATA AB/'SPDFGHIKLMNOQRTUVWXYZ'/

      i = INDEX (AS,a)
      if(i.eq.0) i = INDEX (AB,a)
      if(i.eq.0) i = ICHAR(a)-101
      LA = i-1

      End Function LA


!====================================================================
      Character(4) Function ELF4(n,l,k)
!====================================================================
!     gives the A4 specroscopic notation for electron orbital (n,l,k)
!
!     n must be < 100; if they > 100, we choose 'k' for n=107
!     or 'n' otherwise
!
!     k must be < 61*61; if k<=61 - one character from ASET
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,l,k
      Integer :: i,k1,k2
      Character(4) :: EL
      Character(1), external :: AL

      Character(61) :: ASET = &
       '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

      EL='    '; i=4

      ! index value:

      if(k.gt.0) then
       if(k.le.61) then
        EL(i:i)=ASET(k:k); i=i-1
       elseif(k.le.61*61) then
        k1=k/61; k2=mod(k,61);
        if(k2.eq.0) then; k1=k1-1; k2=61; end if
        EL(i:i)=ASET(k2:k2); i=i-1
        EL(i:i)=ASET(k1:k1); i=i-1
       else
        write(*,*) ' ELF4: set index is out of limits:',k
        Stop
       end if
      end if

      ! l-value

      EL(i:i)=AL(l,1);  i=i-1

      ! n-value

      if(n.gt.0.and.n.lt.10) then
       write(EL(i:i),'(i1)') n
      elseif(n.ge.10.and.n.lt.100.and.i.ge.2) then
       write(EL(i-1:i),'(i2)') n
      elseif(n.eq.107) then
       EL(i:i) = 'k'   ! ICHAR('k')=107
      else
       EL(i:i) = 'n'   ! ICHAR('n')=110
      end if

      ELF4=EL

      End Function ELF4


!====================================================================
      Character(3) Function ELF3(n,l,k)
!====================================================================
!     gives the A3 specroscopic notation for electron orbital (n,l,k)
!
!     n and k must be < 10; if they > 9, they are replaced
!     by characters with code i+ICHAR(1)-1.
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n,l,k
      Character(1), external :: AL

      Character(61) :: ASET = &
       '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

      if(k.eq.0) then
       if(n.lt.100) then
        write(ELF3(1:2),'(i2)') n
       else
        ELF3(1:1) = ' '
        ELF3(2:2)=CHAR(n)
       end if
       ELF3(3:3)=AL(l,1)
      else
       if(n.lt.10) then
        write(ELF3(1:1),'(i1)') n
       else
        ELF3(1:1)=CHAR(n)
       end if
       ELF3(2:2)=AL(l,1)
       ELF3(3:3)=ASET(k:k)
      end if

      End Function ELF3


!====================================================================
      Subroutine EL4_NLK(EL,n,l,k)
!====================================================================
!     decodes the A4 specroscopic notation for orbital (n,l,k)
!
!     It is allowed following notations: 1s, 2s3, 2p30, 20p3,
!     1sh, 20sh, kp, kp1, kp11, ns, ns3, ns33, ... .
!
!     Call:  LA, ICHAR
!--------------------------------------------------------------------
      Implicit none
      Character(4), intent(in) :: EL
      Integer, intent(out) :: n,l,k
      Integer :: i,j,ic, k1,k2
      Integer, external :: LA

      Character(61) :: ASET = &
       '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

      n=0; l=-1; k=0; i=1; j=1
    1 if(EL(i:i).eq.' ') then
       i=i+1
       if(i.le.4) go to 1
      end if
      if(i.eq.4.and.j.eq.1) j=2

      if(i.gt.4) Return

      ic=ICHAR(EL(i:i))-ICHAR('1')+1

      if(j.eq.1) then                      !  n -> ?

       if(n.eq.0.and.ic.le.9) then
        n=ic; j=1
       elseif(n.eq.0.and.ic.gt.9) then
        n=ICHAR(EL(i:i)); j=2
       elseif(n.gt.0.and.ic.le.9) then
        n=n*10+ic; j=2
       else
        l=LA(EL(i:i)); j=3
       end if

      elseif(j.eq.2) then                   !  l -> ?
       l=LA(EL(i:i)); j=3

      elseif(j.eq.3) then                   !  k -> ?
       if(i.eq.3) then
        k1 = INDEX(ASET,EL(i:i)); k2 = INDEX(ASET,EL(i+1:i+1))
        if(k2.gt.0) k = k1*61+k2
        if(k2.eq.0) k = k1
       else
        k = INDEX(ASET,EL(i:i))
       end if
       j=4
      end if

      i=i+1
      if(i.le.4.and.j.le.3) go to 1
      if(n.eq.0.or.l.lt.0) then
       write(*,*) 'EL4_nlk is fail to decode: ',EL
       Stop ' '
      end if

      End Subroutine EL4_NLK


!====================================================================
      Subroutine EL3_NLK(EL,n,l,k)
!====================================================================
!     decodes the A3 specroscopic notation for orbital (n,l,k)
!
!     It is allowed following notations: 1s, 2s3, 20p3, 1sh, kp, kp1,
!                                        ns, ns3, nsa, ... .
!     Call:  LA, ICHAR
!--------------------------------------------------------------------
      Implicit none
      Character(3), intent(in) :: EL
      Integer, intent(out) :: n,l,k
      Integer :: i,j,ic
      Integer, external :: LA

      Character(61) :: ASET = &
       '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

      n=0; l=-1; k=0; i=1
      j=1
    1 if(EL(i:i).eq.' ') then
       i=i+1
       if(i.le.3) go to 1
      end if
      if(i.gt.3) Return

      ic=ICHAR(EL(i:i))-ICHAR('1')+1

      if(j.eq.1) then            ! search n

       if(n.eq.0.and.ic.le.9) then
        n=ic; j=1
       elseif(n.eq.0.and.ic.gt.9) then
        n=ICHAR(EL(i:i)); j=2
       elseif(n.gt.0.and.ic.le.9) then
        n=n*10+ic; j=2
       else
        l=LA(EL(i:i)); j=3
       end if

      elseif(j.eq.2) then        ! search l
       l=LA(EL(i:i));  j=3

      elseif(j.eq.3) then        ! search k

       k = INDEX(ASET,EL(i:i))

      end if

      i=i+1
      if(i.le.3.and.j.le.3) go to 1

      if(n.ge.100.or.l.lt.0.or.k.ge.62) &
       write(*,*) 'EL3_NLK is fail to decode ',EL

      End Subroutine EL3_NLK



