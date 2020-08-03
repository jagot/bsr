!======================================================================
      Subroutine Gen_jj_states (AF_inp,AF_out,jmin,jmax)
!======================================================================
!     preparation the ASF list from a list of (nlj)^q configurations
!----------------------------------------------------------------------
      Use conf_jj; Use orb_jj

      Character(*), intent(in) :: AF_inp, AF_out
      Integer :: jmin, jmax, nu1,nu2, i,ii,ic
      Integer, external :: Icheck_file

! ... Check the files:

      if(Icheck_file(AF_inp).eq.0) Return 
      Call Find_free_unit(nu1)
      open(nu1,file=AF_inp); rewind(nu1)
      Call Find_free_unit(nu2)
      open(nu2,file=AF_out); rewind(nu2)

! ... reload the header:

      Do
       read(nu1,'(a)',end=10) AS; write(nu2,'(a)') trim(AS)
       if(AS.eq.'CSF(s):') Exit
      End do

! ... check jmin and jmax:

      if(jmin.eq.-1.or.jmax.eq.-1) then
       Do 
        read(nu1,'(a)',end=10) AS
        if(AS(6:6).eq.'(') Exit
       End do
       CONFIG = AS;  Call Decode_confj
       i = SUM(iq(1:no)*jn(1:no))
      end if
      if(jmin.eq.-1) jmin=mod(i,2)
      if(jmax.eq.-1) jmax = i
      if(jmax.lt.jmin) jmax=jmin

! ... read configurations for each J-total:

      Do j = jmin,jmax,2;   J_min=j; J_max=j

      ncfg=0
      rewind(nu1)
    1 read(nu1,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(6:6).ne.'(') go to 1

      CONFIG = AS;  Call Decode_confj

      Call Sum_Term

      go to 1
    2 Rewind(nu1)

      if(ncfg.eq.0) Cycle

      Do ic=1,ncfg
       Call Get_cfg_jj(ic)
       Call Incode_cj
       ii = no*9
       write(nu2,'(a)') CONFIG(1:ii)
       write(nu2,'(a)') SHELLJ(1:ii)
       write(nu2,'(9x,a)') INTRAJ(1:ii)
      End do
      if(j.lt.jmax) write(nu2,'(a)') ' *'
      if(j.eq.jmax) write(nu2,'(a)') '* '

      End do ! over j

      BACKSPACE(nu2);  write(nu2,'(a)') '* '

   10 Return

      End  Subroutine Gen_jj_states


!----------------------------------------------------------------------
      Subroutine Sum_Term
!----------------------------------------------------------------------
!     exhaustion of shell-terms
!----------------------------------------------------------------------
      Use conf_jj; Use orb_jj

      Implicit none
      Integer :: mt(msh),nt(msh)
      Integer :: i,ii,i1,i2, JT,JV,JW,JQ
      Integer, external :: Jterm

!     mt(i) - the number of term in shell i
!     nt(i) - the term inder consideration

      i1=1                     ! i1 - low  limit of shells
      i2=no                    ! i2 - high limit of shells
      Do i=i1,i2
       mt(i)=Jterm(jn(i),iq(i),-1,JT,JV,JW,JQ)
      End do

      i=i1                     ! first shell under consideration
      nt = 1

    1 Continue

      ii = Jterm(jn(i),iq(i),nt(i),Jshell(i),Vshell(i),JW,JQ)

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

    3 Continue

      End Subroutine Sum_Term


!----------------------------------------------------------------------
      Subroutine Sum_Iterm
!----------------------------------------------------------------------
!     exhaustion of intermediate terms
!----------------------------------------------------------------------
      USE conf_jj; Use orb_jj
 
      Integer :: js_min(msh),js_max(msh)

      Jintra(1)=Jshell(1)
      if(no.eq.1) then
       if(Jshell(no).ge.J_min.and.Jshell(no).le.J_max)  ic=Iadd_cfg_jj('detect')
       Return
      end if

      i1=2                         ! i1 - low  limit
      i2=no                        ! i2 - high limit in array LS(...)

      i=i1
    1 j1=i-1; j2=i

      js_min(i)=IABS(Jintra(j1)-Jshell(j2))
      js_max(i)=     Jintra(j1)+Jshell(j2) 
      Jintra(i)=js_min(i)

    2 if(i.lt.i2) then
       i=i+1; go to 1
      else

       if(Jintra(no).ge.J_min.and.Jintra(no).le.J_max)  ic=Iadd_cfg_jj('detect')

      end if

    3 if(Jintra(i).lt.js_max(i)) then
       Jintra(i)=Jintra(i)+2
       go to 2
      else
       if(i.eq.i1) go to 4
       i=i-1; go to 3
      end if

    4 Continue
       
      End  Subroutine Sum_Iterm




