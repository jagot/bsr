!======================================================================
!      target_jj,  zf_res   -->  ak_tab
!
!      creates the list of branching ratios (file 'ak_tab') needed to
!      calculate the cascade contribution
!
!======================================================================
       Use target_jj

       IMPLICIT REAL(8) (A-H,O-Z)

       Character(80) :: AS
       Character(13) :: lab,lab1,lab2

       Character(20), Allocatable, Dimension(:) :: LABEL, LABEL1,LABEL2
       REAL(8), Allocatable, DIMENSION(:) :: ET, ET1,ET2
       REAL(8), Allocatable, DIMENSION(:) :: FL,FV,SL,SV,ATL,ATV,AL,AV
       INTEGER, Allocatable, DIMENSION(:) :: jot, jot1,jot2
       INTEGER, Allocatable, DIMENSION(:) :: IPT,IPT1,IPT2

       Real(8), parameter :: Eps_E = 2.D-6
       Logical :: EX
       Character(40) :: AF_tar='target_jj'
       Character(40) :: AF_zf ='zf_res'

       Integer(4) :: mt = 20000

      iarg = IARGC()
      if(iarg.gt.0) then
       write(*,*) 'ak_tab_jj creates the list of branching ratios '
       write(*,*) 'needed to calculate the cascade contribution:'
       write(*,*) '      target_jj,  zf_res   -->  ak_tab '
       Stop ' '
      end if

!----------------------------------------------------------------------
! ... target information:

      Call Check_file(AF_tar)

      nut=1; Open(nut,file=AF_tar);
      Call Read_target_jj(nut)
      close(nut)

      Allocate(ATL(ntarg),ATV(ntarg),LABEL(ntarg))

      Z = nz;  AWT = 0.d0
      Call Conv_au (Z,AWT,au_cm,au_eV,0)
      Ry = au_eV/2.d0

!----------------------------------------------------------------------
! ... zf - data ...

       Call Check_file(AF_zf)
       in=1; Open(in,file=AF_zf,status='OLD')


       Allocate ( jot1(mt),jot2(mt), &
                  IPT(mt),IPT1(mt), IPT2(mt), &
                  FL(mt),FV(mt), SL(mt),SV(mt), AL(mt),AV(mt),&
                  ET1(mt),ET2(mt), LABEL1(mt),LABEL2(mt))

       rewind(in)
       n=0
     3 read(in,'(a)',end=4) AS; if(AS(10:10).ne.'.') go to 3
       read(AS,'(i2,f16.8,2x,a)') j1,E1,LAB1
       read(in,'(i2,f16.8,2x,a)') j2,E2,LAB2
       read(in,'(a)') AS
       read(in,'(8x,D13.5,7x,D13.5,8x,D13.5)') xsl,xfl,xal
       read(in,'(8x,D13.5,7x,D13.5,8x,D13.5)') xsv,xfv,xav

        i1 = 0
        Do i = 1,ntarg
         if(abs(E1-Etarg(i)).gt.Eps_E) Cycle
         if(j1.ne.jtarg(i)) Cycle
         i1 = i
         Exit
        End do
        i2 = 0
        Do i = 1,ntarg
         if(abs(E2-Etarg(i)).gt.Eps_E) Cycle
         if(j2.ne.jtarg(i)) Cycle
         i2 = i
         Exit
        End do


        if(i1.eq.0) go to 3
        if(i2.eq.0) go to 3

        n = n +  1; if(n.gt.mt) Stop ' nt > mt !'

        IPT1(n) = i1; IPT2(n) = i2
        jot1(n) = j1; jot2(n) = j2
        FL(n)   = xfl;  FV(n) = xfv
        SL(n)   = xsl;  SV(n) = xsv
        AL(n)   = xal;  AV(n) = xav
        ET1(n)  = E1;  ET2(n) = E2
        LABEL1(n) = LAB1;  LABEL2(n) = LAB2
        LABEL(i1) = LAB1;  LABEL(i2) = LAB2

        go to 3
    4 Continue

      nt = n

!----------------------------------------------------------------------
! ...  output total AK, life-times and branching ratios in ak_tab:

       ATL = 0.d0; ATV = 0.d0
       Do n = 1,nt
        i=IPT2(n)
        ATL(i) = ATL(i) + AL(n)
        ATV(i) = ATV(i) + AV(n)
       End do

       iout=2; Open(iout,file='ak_tab')

       Call sortI(nt,IPT1,IPT)

       write(iout,'(/a,i5/)') 'nt = ',nt
       Do j=1,nt; i=IPT(j)
        i1=IPT1(i); i2=IPT2(i)
        S = abs(FL(i)-FV(i))/(FL(i)+FV(i))*200
        write(iout,'(2i6,2f10.5,3x,a12,a3,a12,1P2E12.4,0Pf8.1)') &
                  i1,i2, AL(i)/ATL(i2), AV(i)/ATV(i2), &
                  LABEL1(i),' - ',LABEL2(i), FL(i),FV(i),S
       End do

       write(iout,'(/a,i5/)') 'ntarg = ',ntarg
       Do i=1,ntarg
        write(iout,'(a20,5x,a20,f16.8,F10.3,4D13.5)') AFT(i),BFT(i), &
                     Etarg(i), (Etarg(i)-Etarg(1))*au_eV,  &
                     ATL(i),ATV(i), 1.d+9/ATL(i), 1.d+9/ATV(i)
       End do

      Close(iout)

      End ! utility ak_tab

