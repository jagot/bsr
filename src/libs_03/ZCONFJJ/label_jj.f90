!======================================================================
      Subroutine label_jj (mlab,Lab,iset)
!======================================================================
!     Lab  ->   packed configuration notation
!     iset = 0 -> all set numbers forced to 0
!----------------------------------------------------------------------
      Use conf_jj,  mmm => mlab

      Implicit none
      Integer, intent(in) :: mlab,iset
      Character(mlab), intent(out) ::  Lab
      Character(5) :: EL
      Character(5), external :: ELi
      Integer :: i,ii,j,jj,k,kk, JT,JV,JW,JQ
      Integer, external :: Jterm, Jparity

      Lab = ' '

      k=1
      Do i=1,no; if(iq(i).lt.1) Cycle

       kk = jn(i)+1
       if(k.eq.1.and.no-i.ge.3.and.iq(i).eq.kk) Cycle

       if(i.eq.1.and.no.ge.2.and.iq(i).eq.kk.and.ln(i).eq.0) Cycle


       ! ... orbital

       ii=in(i); if(iset.eq.0) ii=0      
       EL=ELi(nn(i),kn(i),ii)
       Do j = 1,5
        if(EL(j:j).eq.' ') Cycle
        Lab(k:k)=EL(j:j); k=k+1
       End do

       ! ... number of electrons

       if(iq(i).gt.1.and.iq(i).le.9) then
        write(Lab(k:k),'(i1)') iq(i); k=k+1
       elseif(iq(i).gt.9) then
        write(Lab(k:k+1),'(i2)') iq(i); k=k+2
       end if
       Lab(k:k)='.'; k=k+1

       ! ... shell term
      
       if(Jterm (jn(i),iq(i),-1,JT,JV,JW,JQ).gt.1) then

        k=k-1;  Lab(k:k)='_'; k=k+1

        kk = mod(Jshell(i),2)
        if(kk.eq.0) then
         jj = Jshell(i)/2
         if(jj.le.9) then
          write(Lab(k:k),'(i1)') jj; k=k+1
         else
          write(Lab(k:k+1),'(i2)') jj; k=k+2
         end if
        else
         jj = Jshell(i)
         if(jj.le.9) then
          write(Lab(k:k),'(i1)') jj; k=k+1
         else
          write(Lab(k:k+1),'(i2)') jj; k=k+2
         end if
          write(Lab(k:k+1),'(a2)') '/2'; k=k+2
        end if

        Lab(k:k)='.'; k=k+1

       end if

       ! ... intermediate term

       if(i.eq.1.and.no.gt.1) Cycle
       kk = 0 
       if(i.gt.1) then
        if(iabs(Jshell(i)-Jintra(i-1)).lt.Jshell(i)+Jintra(i-1)) kk=1
       end if

       if(kk.eq.1.or.i.eq.no) then
!          i.eq.no) then
        kk = mod(Jintra(i),2)
        if(kk.eq.0) then
         jj = Jintra(i)/2
         if(jj.le.9) then
          write(Lab(k:k),'(i1)') jj; k=k+1
         else
          write(Lab(k:k+1),'(i2)') jj; k=k+2
         end if
        else
         jj = Jintra(i)
         if(jj.le.9) then
          write(Lab(k:k),'(i1)') jj; k=k+1
         else
          write(Lab(k:k+1),'(i2)') jj; k=k+2
         end if
          write(Lab(k:k+1),'(a2)') '/2'; k=k+2
        end if

        if(i.ne.no) then; Lab(k:k)='_'; k=k+1; end if
       end if

      if(k-1.gt.mlab) write(*,*) 'Label_jj too long'

      End do

      i = Jparity(no,ln,iq); if(i.eq.-1) LAB(k:k)='*'

      End Subroutine label_jj


!======================================================================
      Subroutine R_label_jj(nuc,kset)
!======================================================================
!     generates configuration LABEL list from c-file (unit nuc)
!     LABEL list is located in the module conf_jj   
!     kset=0 -> nulify all set numbers
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer, intent(in) :: nuc,kset
      Integer ::  nc

      if(ncfg.eq.0) Return
      if(allocated(LABEL)) Deallocate(LABEL);  Allocate(LABEL(ncfg))
      nc=0
      rewind(nuc)
    1 read(nuc,'(a)',end=2) CONFIG
      if(CONFIG(1:3).eq.'***') go to 2
      if(CONFIG(6:6).ne.'(') go to 1
      read(nuc,'(a)') SHELLJ
      read(nuc,'(5x,a)') INTRAJ
      Call Decode_cj
      nc=nc+1
      if(nc.gt.ncfg) Stop 'R_label: nc  > ncfg'
      Call Label_jj (mlab,Label(nc),kset)
      go to 1
    2 Continue
      if(nc.ne.ncfg) Stop ' R_label: nc <> ncfg'

      End Subroutine R_label_jj



