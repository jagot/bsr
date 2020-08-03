!======================================================================
      Subroutine Incode_c
!======================================================================
!     incodes the configuration from ZOI format to c-file format
!----------------------------------------------------------------------
      Use conf_LS

      Implicit none
      
      Character(4), External :: ELF4
      Character(1), External :: AL
      Character(4) :: BLANK = '    '
      Integer(4) :: i,k,m

      k=0; m=1
      Do i=1,no; k=k+iq(i); if(iq(i).eq.0) m=0; End do
      
      if(m.eq.0.and.k.gt.0) then
       k=0
       Do i=1,no
        if(iq(i).eq.0) Cycle 
        k=k+1
        nn(k)=nn(i); ln(k)=ln(i); kn(k)=kn(i); iq(k)=iq(i)
        LS(k,:)=LS(i,:)
       End do         
       no=k
      end if

       k=0
       Do i=1,no
        CONFIG(k+1:k+4) = ELF4(nn(i),ln(i),kn(i))
        write(CONFIG(k+5:k+8),'(a1,i2,a1)') '(',iq(i),')'
        k=k+8
       End do

       Do i=no+1,msh
        CONFIG(k+1:k+4)=BLANK
        CONFIG(k+5:k+8)=BLANK
        k=k+8
       End do

       k=0
       Do i=1,no
        write(COUPLE(k+1:k+4),'(i2,a1,i1)')  &
                              LS(i,3),AL(LS(i,2),6),LS(i,1)
        k=k+4
       End do

       Do i=2,no
        write(COUPLE(k+1:k+4),'(i2,a1,a1)')  &
                              LS(i,5),AL(LS(i,4),6),' '
        k=k+4
       End do

       Do i=no+no,msh+msh
        COUPLE(k+1:k+4)=BLANK
        k=k+4
       End do

       End Subroutine Incode_c


!======================================================================
       Subroutine Decode_c
!======================================================================
!      decode the config. from c-file format to ZOI format
!----------------------------------------------------------------------

       USE conf_LS

       no=0
       ii = LEN_TRIM(CONFIG)
       ii = ii/8
       k=1
       j=2
       Do i=1,ii
        if(CONFIG(k+4:k+4).ne.'(') Exit
        no=no+1
        Call EL4_nlk(CONFIG(k:k+3),nn(i),ln(i),kn(i))
        read(CONFIG(k+5:k+6),'(i2)') iq(i)
        k=k+8
        read(COUPLE(j:j),'(i1)') LS(i,3)
        LS(i,2)=2*LA(COUPLE(j+1:j+1))+1
        read(COUPLE(j+2:j+2),'(i1)') LS(i,1)
        j=j+4
       End do

       LS(1,4)=LS(1,2)
       LS(1,5)=LS(1,3)
       Do i=2,no
        read(COUPLE(j:j),'(i1)') LS(i,5)
        LS(i,4)=2*LA(COUPLE(j+1:j+1))+1
        j=j+4
       End do

       Ltotal = LS(no,4)
       Stotal = LS(no,5)

       End Subroutine Decode_c

!======================================================================
       Subroutine Decode_cn(CONFIG,COUPLE,no,nn,ln,iq,kn,LS)
!======================================================================
!      decode the config. from c-file format to ZOI format
!----------------------------------------------------------------------
       Implicit none
       Character(*), intent(in) :: CONFIG,COUPLE
       Integer :: no,nn(*),ln(*),iq(*),kn(*),LS(5,*)
       Integer :: ii, k,i,j 
       Integer, external :: LA

       no=0; ii=LEN_TRIM(CONFIG); ii=ii/8;  k=1; j=2
       Do i=1,ii
        if(CONFIG(k+4:k+4).ne.'(') Exit
        no=no+1
        Call EL4_nlk(CONFIG(k:k+3),nn(i),ln(i),kn(i))
        read(CONFIG(k+5:k+6),'(i2)') iq(i)
        k=k+8
        read(COUPLE(j:j),'(i1)') LS(3,i)
        LS(2,i)=2*LA(COUPLE(j+1:j+1))+1
        read(COUPLE(j+2:j+2),'(i1)') LS(1,i)
        j=j+4
       End do

       LS(4,1)=LS(2,1)
       LS(5,1)=LS(3,1)
       Do i=2,no
        read(COUPLE(j:j),'(i1)') LS(5,i)
        LS(4,i)=2*LA(COUPLE(j+1:j+1))+1
        j=j+4
       End do

       End Subroutine Decode_cn


!====================================================================
      Subroutine Incode_conf(ASYM,BSYM)
!====================================================================
!
!     incodes the conf. and angular symmetries (ASYM and BSYM)
!     (ASYM also includes the total term!)
!--------------------------------------------------------------------

      Use conf_LS
      IMPLICIT NONE

      Character(26), Intent(out) :: ASYM
      Character(40), Intent(out) :: BSYM
      Character(1), EXTERNAL :: AL
 
      Integer(4) :: i,k,m

      k=1
      m=1
      Do i=1,no
       write(ASYM(k:k+2),'(a1,i2)') AL(ln(i),1),iq(i)
       k=k+3
       write(BSYM(m:m+4),'(2a1,i1,2a1)') &
                    AL(LS(i,3),2),AL(LS(i,2),6),LS(i,1), &
                    AL(LS(i,5),2),AL(LS(i,4),6)
       m=m+5
      End do

      Do i=no+1,msh
       ASYM(k:k+2)='   '
       k=k+3
       BSYM(m:m+4)='     '
       m=m+5
      End do

      ASYM(25:25) = AL(LS(no,4),6)
      ASYM(26:26) = AL(LS(no,5),2)

      End Subroutine Incode_conf



!====================================================================
      Subroutine Decode_ASYM(ASYM)
!====================================================================
!     decode conf. symmetry ASYM which also includes the total term
!--------------------------------------------------------------------

      Use conf_LS
      Character(26), Intent(in) :: ASYM
 
      no = LEN_TRIM(ASYM(1:24))/3
      k=1
      Do i=1,no
       ln(i)=LA(ASYM(k:k))
       read(ASYM(k+1:k+2),*) iq(i)
       k=k+3
      End do
      Ltotal = LA(ASYM(25:25)); Ltotal=Ltotal+Ltotal+1
      Stotal = LA(ASYM(26:26))

      End Subroutine Decode_ASYM


!====================================================================
      Subroutine Decode_BSYM(BSYM)
!====================================================================
!     decodes the angular symmetry BSYM
!--------------------------------------------------------------------

      Use conf_LS

      Character(40), Intent(in) :: BSYM

      no = LEN_TRIM(BSYM)/5
      m=1
      Do i=1,no
       LS(i,3) = LA(BSYM(m:m))
       LS(i,2) = LA(BSYM(m+1:m+1)); LS(i,2)=LS(i,2)+LS(i,2)+1
       read(BSYM(m+2:m+2),*) LS(i,1)
       LS(i,5) = LA(BSYM(m+3:m+3))
       LS(i,4) = LA(BSYM(m+4:m+4)); LS(i,4)=LS(i,4)+LS(i,4)+1
       m=m+5
      End do

      End Subroutine Decode_BSYM


!======================================================================
      Subroutine Conf_string(no,nn,ln,iq,label)
!======================================================================
!     configuration label from ZOI format to c-file format
!----------------------------------------------------------------------

      Implicit none
      
      Integer, Intent(in) :: no,nn(no),ln(no),iq(no)
      Character(*) :: label

      Character(4), External :: ELF4
      Character :: BLANK = ' '
      Integer :: i,k,m

      label = BLANK

      k=0
      Do i=1,no; if(iq(i).eq.0) Cycle
       label(k+1:k+4) = ELF4(nn(i),ln(i),0)
       write(label(k+5:k+8),'(a1,i2,a1)') '(',iq(i),')'
       k=k+8
      End do

      Call Clean_a(label)

      End Subroutine Conf_string



!======================================================================
      Subroutine Conf_label(no,nn,ln,iq,label,kk)
!======================================================================
!     configuration label from ZOI format to c-file format
!----------------------------------------------------------------------

      Implicit none
      
      Integer, Intent(in) :: no,nn(no),ln(no),iq(no),kk
      Character(*) :: label

      Character(4), External :: ELF4
      Character :: BLANK = ' '
      Integer(4) :: i,k,m

      label = BLANK

      m=0; if(no.le.kk) m=1
      k=0
      Do i=1,no; if(iq(i).eq.0) Cycle
       if(i.eq.no) m = 1
       if(m.eq.0.and.i.gt.no-kk) m = 1
       if(m.eq.0.and.4*ln(i)+2.eq.iq(i)) Cycle
       m = 1
       label(k+1:k+4) = ELF4(nn(i),ln(i),0)
       k=k+4
       if(iq(i).gt.1) then
        write(label(k+1:k+2),'(i2)') iq(i)
        k=k+2
       end if
       k=k+1; write(label(k:k),'(a)') '_'
      End do
      
      label(k:k)=BLANK

      Call Clean_a(label)

      End Subroutine Conf_label


!======================================================================
       Subroutine Decode_conf
!======================================================================
!      decode the config. from c-file format to ZOI format
!----------------------------------------------------------------------
       Use conf_LS

       no=0
       ii = LEN_TRIM(CONFIG)
       ii = ii/8
       k=1
       Do i=1,ii
        if(CONFIG(k+4:k+4).ne.'(') Exit
        no=no+1
        Call EL4_nlk(CONFIG(k:k+3),nn(i),ln(i),kn(i))
        read(CONFIG(k+5:k+6),'(i2)') iq(i)
        k=k+8
       End do

       End Subroutine Decode_conf
