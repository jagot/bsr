!=====================================================================
!     utility     Z G E N T E R M
!
!                 C O P Y R I G H T -- 2004
!
!     Written by:   Oleg Zatsarinny
!======================================================================
!
!     zgenterm prepares list of states from list of configurations:
!
!       conf.inp --> cfg.inp
!
!     Interactive input/output    or
!
!       zgenconf conf-file J_min J_max L_min L_max S_min S_max c-file
!
!     All ang.momentum values in (2J+1) representation
!     (0 - means 'no restrictions in this respect' )
!
!     Note: there is restriction on max. 2L+1 = 13 for  NONH
!--------------------------------------------------------------------
      Use conf_LS

      Implicit none
      Character(80) :: AS
      Integer :: nu1=1; Character(40) :: AF_inp = 'conf.inp'
      Integer :: nu2=2; Character(40) :: AF_out = 'cfg.inp'
      Integer :: i,i1,i2, iarg

      iarg = command_argument_count()
      if(iarg.gt.0)  Call get_command_argument(1,AF_inp)

      if(AF_inp.eq.'?') then
        write(*,'(/a)') 'zgenterm prepares list of states from list of configurations:'         
        write(*,'(/a)') '         conf.inp --> cfg.inp '                                          
        write(*,'(/a)') 'Interactive input/output    or '                                
        write(*,'(/a)') 'zgenconf conf-file J_min J_max L_min L_max S_min S_max c-file'
        write(*,'(/a)') 'All ang.momentum values in (2J+1) representation'               
        write(*,'(/a)') '(0 - means no restrictions in this respect)'                 
        write(*,'(/a)') 'There is restriction on max. 2L+1 = 13 in accordance with NONH' 
        write(*,'(/a)') 'Example: zgenterm a.conf 0 0 3 3 1 3 a.c' 
        Stop ' '         
      end if             

      if(iarg.ge.8)  Call get_command_argument(8,AF_out)

! ... files

      open(nu1,file=AF_inp,status='OLD')
      open(nu2,file=AF_out)

!------------------------------------------------------------------
! ... input data:

      if(iarg.ge.7) then

       Call get_command_argument(2,AS); read(AS,*) J_min  
       Call get_command_argument(3,AS); read(AS,*) J_max  
       Call get_command_argument(4,AS); read(AS,*) L_min  
       Call get_command_argument(5,AS); read(AS,*) L_max  
       Call get_command_argument(6,AS); read(AS,*) S_min  
       Call get_command_argument(7,AS); read(AS,*) S_max  
 
      else

       write(*,*) ' Enter the J_min,J_max (2J+1):'
       read(*,*) J_min,J_max
       write(*,*) ' Enter the L_min,L_max (2L+1):'
       read(*,*) L_min,L_max
       write(*,*) ' Enter the S_min,S_max (2S+1):'
       read(*,*) S_min,S_max

      end if

!--------------------------------------------------------------------
! ... read the configurations:

      ncfg=0
      Call alloc_cfg_LS(icfg)
      rewind(nu1)
      rewind(nu2)
      read(nu1,'(a)') AS; write(nu2,'(a)') AS
      read(nu1,'(a)') AS; write(nu2,'(a)') AS
    1 read(nu1,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1

      CONFIG = AS(1:64)
      Do i=1,16+16-1
       i1=(i-1)*4+1; i2=i1+3
       COUPLE(i1:i2) = ' 1S0'
      End do
      conf_mode = 0
      if(INDEX(AS,'*').ne.0) conf_mode=1

      Call Decode_c

      Call Sum_Term

      go to 1
    2 Continue

!      Call Test_a

      Do i=1,ncfg
       Call Pri_conf (nu2,i,0.d0) 
      End do
      write(nu2,'(a)') '*'

      write(*,'(a,a,a,a,i6)')  trim(AF_inp),' -> ', trim(AF_out),'  ncfg =',ncfg

      End  ! utility     Z G E N T E R M



!----------------------------------------------------------------------
      Subroutine Sum_Term
!----------------------------------------------------------------------
!     exhaustion of shell-terms
!----------------------------------------------------------------------

      Use conf_LS, only: msh,no,ln,iq,LS

      Implicit none 
      Integer :: mt(msh),nt(msh)
      Integer :: i,i1,i2,ii,IA,IL,IS
      Integer, external :: Iterm_LS

!     mt(i) - the number of term in shell i
!     nt(i) - the term inder consideration

      i1=1                     ! i1 - low  limit of shells
      i2=no                    ! i2 - high limit of shells
      Do i=i1,i2
       mt(i)=Iterm_LS(ln(i),iq(i),-1,IA,IL,IS)
      End do

      i=i1                     ! first shell under consideration
      nt(i)=1

    1 ii=Iterm_LS(ln(i),iq(i),nt(i),LS(i,1),LS(i,2),LS(i,3))
      if(i.lt.i2) then
         i=i+1; nt(i)=1; go to 1
      else
         CALL Sum_Iterm
      end if

    2 nt(i)=nt(i)+1
      if(nt(i).gt.mt(i)) then
        if(i.eq.i1) go to 3
        i=i-1; go to 2
        end if
      go to 1

    3 Return

      End  ! Subroutine Sum_Term


!----------------------------------------------------------------------
      Subroutine Sum_Iterm
!----------------------------------------------------------------------
!     exhaustion of intermediate terms
!----------------------------------------------------------------------
      Use conf_LS, only: msh,no,LS,conf_mode

      Implicit none
      Integer :: LL_min(msh),LL_max(msh),SS_min(msh),SS_max(msh)
      Integer :: i,i1,i2,j1,j2

      LS(1,4)=LS(1,2)
      LS(1,5)=LS(1,3)
      if(conf_mode.eq.1) LS(1,5) = -1
      if(no.eq.1) then;  CALL Output_c; Return; end if

      i1=2                         ! i1 - low  limit
      i2=no                        ! i2 - high limit in array LS(...)
      i=i1
    1 j1=i-1; j2=i

      LL_min(i)=IABS(LS(j1,4)-LS(j2,2))+1
      LL_max(i)=     LS(j1,4)+LS(j2,2) -1
      SS_min(i)=IABS(LS(j1,5)-LS(j2,3))+1
      SS_max(i)=     LS(j1,5)+LS(j2,3) -1
      LS(i,4)=LL_min(i)
      LS(i,5)=SS_min(i)

    2 if(i.lt.i2) then
       i=i+1; go to 1
      else
       CALL Output_c
      end if

    3 if(LS(i,5).lt.SS_max(i)) then
         LS(i,5)=LS(i,5)+2
         go to 2
      elseif(LS(i,4).lt.LL_max(i)) then
         LS(i,4)=LS(i,4)+2
         LS(i,5)=SS_min(i)
         go to 2
      else
         if(i.le.i1) go to 4
         i=i-1; go to 3
      end if

    4 Return

      End  ! Subroutine Sum_Iterm


!======================================================================
      Subroutine Output_c
!======================================================================

      USE conf_LS
       
      Implicit none

      Integer :: i, ILT,IST, j1,j2
      Integer, external :: Ifind_cfg_LS
      

! ... check the total term: 

      ILT=LS(no,4)
      IST=LS(no,5)
      j1=iabs(ILT-IST)+1
      j2=iabs(ILT+IST)-1

      if(L_min.gt.0.and.L_min.gt.ILT) Return
      if(L_max.gt.0.and.L_max.lt.ILT) Return
      if(S_min.gt.0.and.S_min.gt.IST) Return
      if(S_max.gt.0.and.S_max.lt.IST) Return
      if(J_min.gt.0.and.J_min.gt. j2) Return
      if(J_max.gt.0.and.J_max.lt. j1) Return

! NONH restriction:
      Do i=1,no; if(LS(i,2).gt.13) Return; End do


! ... add configuration: 

      if(conf_mode.gt.0) then
       LS = 1;    LS(1,5)=-1;  LS(no,4)=ILT; LS(no,5)=IST
      else
       Call Test_c
      end if

      i = Ifind_cfg_LS()

      End ! Subroutine Output_c


