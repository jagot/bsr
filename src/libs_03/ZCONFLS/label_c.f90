!======================================================================
      Subroutine Label_c (Labelc,n1,kset)
!======================================================================
!     Labelc  ->   packed configuration notation
!     n1 -> begin from shell n1
!     kset = 0 -> all set numbers --> 0
!----------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer, intent(in) :: n1,kset
      Character(*) ::  Labelc
      Character(4), external :: ELF4
      Character(1), external :: AL
      Character(4) :: EL

      Integer :: i,i1,i2,j,k,kk

      Labelc = ' '

      k=1; i1=1; if(n1.gt.1.and.n1.le.no) i1=n1

      i2=i1
      Do i=i1,no
       if(iq(i).ne.4*ln(i)+2) Exit
       i2 = i+1
      End do
      i1 = min(i2,no)

      Do i=i1,no
       if(iq(i).lt.1) Cycle
                                       ! orbital
       kk=kn(i); if(kset.eq.0) kk=0
       EL=ELF4(nn(i),ln(i),kk)
       j=1; if(EL(1:1).eq.' ') j=2; if(EL(2:2).eq.' ') j=3
       Labelc(k:k+4-j)=EL(j:4)
       k=k+5-j
                                       ! number of electrons
       if(iq(i).gt.1) then
        Labelc(k:k)='('
        k=k+1
        if(iq(i).le.9) then
         write(Labelc(k:k),'(i1)') iq(i)
         k=k+1
        else
         write(Labelc(k:k+1),'(i2)') iq(i)
         k=k+2
        end if
        Labelc(k:k)=')'
        k=k+1
       end if
                                       ! shell term

       if(iq(i).ne.1 .and. iq(i).ne.4*ln(i)+2 .and.  &
          iq(i).ne.4*ln(i)+1) then

        write(Labelc(k:k+2),'(i1,a1,i1)')            &
                   LS(i,3),AL(LS(i,2),6),LS(i,1)
        k=k+3
        if(LS(i,1).eq.0) k=k-1
       end if
                                       ! intermediate term
       if(i.gt.i1.and.iq(i).ne.4*ln(i)+2) then
        write(Labelc(k:k+2),'(a1,i1,a1)')  &
                   '_',LS(i,5),AL(LS(i,4),6)
        k=k+3
       end if

       if(i.ne.no) write(Labelc(k:k),'(a1)') '_'
       k=k+1
      End do

      End Subroutine Label_c


!======================================================================
      Subroutine Label_cc (Labelc,n1,kset)
!======================================================================
!     Labelc  ->   packed configuration notation
!     n1 -> begin from shell n1
!     kset = 0 -> all set numbers --> 0
!----------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Integer, intent(in) :: n1,kset
      Character(*) ::  Labelc
      Character(4), external :: ELF4
      Character(1), external :: AL
      Character(4) :: EL

      Integer :: i,i1,i2,j,k,kk

      Labelc = ' '

      k=1; i1=1; if(n1.gt.1.and.n1.le.no) i1=n1

      i2=i1
      Do i=i1,no
       if(iq(i).ne.4*ln(i)+2) Exit
       if(no-i.lt.2) Exit
       i2 = i+1
      End do
      i1 = min(i2,no)

      Do i=i1,no
       if(iq(i).lt.1) Cycle
       if(i.lt.no.and.iq(i).eq.4*ln(i)+2) Cycle
                                       ! orbital
       kk=kn(i); if(kset.eq.0) kk=0
       EL=ELF4(nn(i),ln(i),kk)
       j=1; if(EL(1:1).eq.' ') j=2; if(EL(2:2).eq.' ') j=3
       Labelc(k:k+4-j)=EL(j:4)
       k=k+5-j
                                       ! number of electrons
       if(iq(i).gt.1) then
!        Labelc(k:k)='(';   k=k+1
        if(iq(i).le.9) then
         write(Labelc(k:k),'(i1)') iq(i)
         k=k+1
        else
         write(Labelc(k:k+1),'(i2)') iq(i)
         k=k+2
        end if
        Labelc(k:k)='_';    k=k+1
       end if
                                       ! shell term

       if(iq(i).ne.1 .and. iq(i).ne.4*ln(i)+2 .and.  &
          iq(i).ne.4*ln(i)+1) then

        write(Labelc(k:k+2),'(i1,a1,i1)')            &
                   LS(i,3),AL(LS(i,2),6),LS(i,1)
        k=k+3
        !        if(LS(i,1).eq.0) k=k-1
       end if
                                       ! intermediate term
       if(i.gt.i1.and.iq(i).ne.4*ln(i)+2) then
        write(Labelc(k:k+2),'(a1,i1,a1)')  &
                   '_',LS(i,5),AL(LS(i,4),6)
        k=k+3
       end if

       if(i.ne.no) write(Labelc(k:k),'(a1)') '_'
       k=k+1
      End do

      End Subroutine Label_cc
