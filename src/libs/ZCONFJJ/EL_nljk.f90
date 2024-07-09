!=======================================================================
      Subroutine EL_nljk(EL,n,kappa,l,j,k)
!=======================================================================
!     decodes the spectroscopic notation for electron orbital (n,l,j,k)
!
!     It is allowed following notations: 1s , 2s 3, 2p-30, 20p-3,
!     1s h, 20s h, kp , kp 1, kp-11, ns , ns 3, ns 33, ... .
!
!     Call:  LA, INDEX
!----------------------------------------------------------------------
      Use orb_jj, only: kset, ASET

      Implicit none
      Character(5), intent(in) :: EL
      Integer, intent(out) :: n,l,j,k,kappa
      Integer :: jj, k1,k2, n1,n2
      Integer, external :: LA, kappa_lj

      jj=0
      Do j=5,3,-1
       if(EL(j:j).eq.'-'.or.EL(j:j).eq.'+') then; jj=j; Exit; end if
      End do
      if(jj.eq.0) then
       Do j=5,3,-1
        if(EL(j:j).eq.' '.and.EL(j-1:j-1).ne.' ') then
         jj=j; Exit
        end if
       End do
      end if

      if(jj.eq.5) then

       k = 0
       l = LA(EL(4:4))
       n1 = INDEX(ASET,EL(3:3))
       n2 = INDEX(ASET,EL(2:2))

      elseif(jj.eq.4) then

       k = INDEX(ASET,EL(5:5))
       l = LA(EL(3:3))
       n1 = INDEX(ASET,EL(2:2))
       n2 = INDEX(ASET,EL(1:1))

      elseif(jj.eq.3) then

       k1 = INDEX(ASET,EL(5:5))
       k2 = INDEX(ASET,EL(4:4))
       k = k2*kset+k1
       l = LA(EL(2:2))
       n1 = INDEX(ASET,EL(1:1))
       n2 = 0

      else

       write(*,*) 'EL_NLJK: can not decode ',EL
       Stop ' '

      end if

       if(n1.le.9.and.n2.le.9) then
        n = n2*10+n1
       elseif(n2.eq.0.and.n1.eq.20) then
        n=ICHAR('k')
       elseif(n2.eq.0.and.n1.eq.23) then
        n=ICHAR('n')
       else
        n = n2*kset+n1
       end if


       j = l+l+1; if(EL(jj:jj).eq.'-') j = l+l-1
       kappa = kappa_lj(l,j)

      End Subroutine EL_NLJK


!=======================================================================
      Character(5) Function ELi(n,kappa,k)
!=======================================================================
!     provides the spectroscopic notation for electron orbital (n,l,j,k)
!     set index k must be < 61*61 (see ASET for incoding k)
!-----------------------------------------------------------------------
      Use orb_jj, only: kset, ASET

      Implicit none
      Integer :: n,l,j,k, ll,jj, i,k1,k2,n1,n2, kappa
      Character(5) :: EL
      Character(1), external :: AL
      Integer, external :: l_kappa, j_kappa

      l = l_kappa(kappa)
      j = j_kappa(kappa)

      if(n.le.0.or.l.lt.0.or.j.le.0.or.k.lt.0) then
       write(*,'(a,4i5)') &
           ' ELj: parmeters are out of limits: n,l,j,k=',n,l,j,k
       Error Stop 'Stop in Elj'
      end if

      EL='    '; i=5

      if(k.lt.0) Stop 'ELi: set index < 0'
      if(k.gt.0) then
       if(k.le.kset) then
        EL(i:i)=ASET(k:k); i=i-1
       else
        k1=k/kset; k2=mod(k,kset);
        if(k2.eq.0) then; k1=k1-1; k2=kset; end if
        if(k1.gt.kset) Stop 'ELi: set index too big'
        EL(i:i)=ASET(k2:k2); i=i-1
        EL(i:i)=ASET(k1:k1); i=i-1
       end if
      end if

      if(j.eq.l+l-1) EL(i:i) = '-'; i=i-1

      EL(i:i)=AL(l,1);  i=i-1

      if(n.lt.0) Stop 'ELi: n < 0'
      if(n.eq.ICHAR('k')) then
        EL(i:i)='k'
      elseif(n.eq.ICHAR('n')) then
        EL(i:i)='n'
      elseif(n.lt.10) then
        write(EL(i:i),'(i1)') n
      elseif(i.ge.2.and.n.lt.100) then
        write(EL(i-1:i),'(i2)') n
      else
        EL(i:i)='k'
      end if

      ELi = EL

      End Function ELi


!=======================================================================
      Character(5) Function ELj(n,l,j,k)
!=======================================================================
!     provides the spectroscopic notation for electron orbital (n,l,j,k)
!     set index k must be < 61*61 (see ASET for incoding k)
!-----------------------------------------------------------------------
      Use orb_jj, only: kset, ASET

      Implicit none
      Integer :: n,l,j,k, ll,jj, i,k1,k2,n1,n2, kappa
      Character(5) :: EL
      Character(1), external :: AL

      if(n.le.0.or.l.lt.0.or.j.le.0.or.k.lt.0) then
       write(*,'(a,4i5)') &
           ' ELj: parmeters are out of limits: n,l,j,k=',n,l,j,k
       Error Stop 'Stop in Elj'
      end if

      EL='    '; i=5

      if(k.lt.0) Stop 'ELi: set index < 0'
      if(k.gt.0) then
       if(k.le.kset) then
        EL(i:i)=ASET(k:k); i=i-1
       else
        k1=k/kset; k2=mod(k,kset);
        if(k2.eq.0) then; k1=k1-1; k2=kset; end if
        if(k1.gt.kset) Stop 'ELi: set index too big'
        EL(i:i)=ASET(k2:k2); i=i-1
        EL(i:i)=ASET(k1:k1); i=i-1
       end if
      end if

      if(j.eq.l+l-1) EL(i:i) = '-'; i=i-1

      EL(i:i)=AL(l,1);  i=i-1

      if(n.lt.0) Stop 'ELi: n < 0'
      if(n.eq.ICHAR('k')) then
        EL(i:i)='k'
      elseif(n.eq.ICHAR('n')) then
        EL(i:i)='n'
      elseif(n.lt.10) then
        write(EL(i:i),'(i1)') n
      elseif(i.ge.2.and.n.lt.100) then
        write(EL(i-1:i),'(i2)') n
      else
        EL(i:i)='k'
      end if

      ELj = EL

      End Function ELj



!====================================================================
      Character FUNCTION AL(L,k)
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

      End FUNCTION AL


!====================================================================
      Integer FUNCTION LA(a)
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

      End FUNCTION LA

