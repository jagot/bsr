!======================================================================
      Subroutine Incode_cj
!======================================================================
!     encodes the configuration from INTEGER format to the GRASP
!     c-file format (see module conf_jj for variables)
!----------------------------------------------------------------------
      Use conf_jj
      Implicit none
      Integer :: i,j,m,k,j1,j2,jj
      Character(5), external :: ELi

      CONFIG=' ';  SHELLJ=' ';   INTRAJ=' '

      m=-9; jj=0
      Do i=1,no
       m = m + 9
       CONFIG(m+1:m+5) = ELi(nn(i),kn(i),in(i))
       write(CONFIG(m+6:m+9),'(a1,i2,a1)') '(',iq(i),')'
       SHELLJ(m+1:m+9) = '        '
       INTRAJ(m+1:m+9) = '        '
       if(iq(i).eq.jn(i)+1) Cycle

       k = mod(Jshell(i),2)
       if(k.eq.0) then
        write(SHELLJ(m+1:m+9),'(i9)') Jshell(i)/2
       else
        write(SHELLJ(m+1:m+9),'(i7,a2)') Jshell(i),'/2'
       end if
       if(jn(i).eq.7.and.iq(i).eq.4.and.Jshell(i).gt.0)  &
         write(SHELLJ(m+1:m+5),'(i4,a1)') Vshell(i),';'

       if(i.eq.1.or.i.eq.no) Cycle

       if(SHELLJ(m+1:m+9).ne.'        ') then
        Do j=1,i-1; j1=(j-1)*9+1; j2=j*9
         if(SHELLJ(j1:j2).ne.'        ') jj=1
        End do
       end if

       !if(Jintra(i).eq.0) Cycle
       if(Jshell(i).eq.0) Cycle

       j1 = iabs(Jintra(i-1)-Jshell(i))
       j2 = iabs(Jintra(i-1)+Jshell(i))
       if(j1.eq.j2.and.jj.eq.0) Cycle
       k = mod(Jintra(i),2)
       if(k.eq.0) then
        write(INTRAJ(m+1:m+9),'(i9)') Jintra(i)/2
       else
        write(INTRAJ(m+1:m+9),'(i7,a2)') Jintra(i),'/2'
       end if
       jj = 1
       Call Clean_a(INTRAJ(m+1:m+9))
      End do

      k = mod(Jintra(no),2)
      if(k.eq.0) then
       write(INTRAJ(m+1:m+7),'(i7)') Jintra(no)/2
      else
       write(INTRAJ(m+1:m+7),'(i5,a2)') Jintra(no),'/2'
      end if

      k=0; Do i=1,no; k=k+ln(i)*iq(i); End do; parity=(-1)**k
      if(parity.eq.+1) INTRAJ(m+8:m+9)='+ '
      if(parity.eq.-1) INTRAJ(m+8:m+9)='- '
      Call Clean_a(INTRAJ(m+1:m+9))

      ia = no*9

      End Subroutine Incode_cj


!======================================================================
      Subroutine Decode_cj
!======================================================================
!     decode the configuration from c-file format to INTEGER format
!----------------------------------------------------------------------
      Use conf_jj
      Implicit none
      Integer :: i,j,k,m
      Character(9) :: bl = '        '

      ia=INDEX(CONFIG,')',BACK=.TRUE.); no=ia/9

      Vshell=0; Jshell=0; Jintra=0
      m=-9
      Do i=1,no
       m = m + 9

       Call EL_nljk(CONFIG(m+1:m+5),nn(i),kn(i),ln(i),jn(i),in(i))
       read(CONFIG(m+7:m+8),'(i2)') iq(i)

       if(SHELLJ(m+5:m+5).eq.';') then
        read(SHELLJ(m+4:m+4),'(i1)') Vshell(i)
        read(SHELLJ(m+9:m+9),'(i1)') J; Jshell(i) = 2*J
       elseif(SHELLJ(m+8:m+8).eq.'/') then
        read(SHELLJ(m+1:m+7),'(i7)') Jshell(i)
       elseif(SHELLJ(m+9:m+9).ne.' ') then
        read(SHELLJ(m+1:m+9),'(i9)') J; Jshell(i) = 2*J
       end if

       k = INDEX(INTRAJ(m+1:m+9),'/')
       if(i.eq.1) then
        Jintra(i) = Jshell(i)
       elseif(i.eq.no) then
        Cycle
       elseif(k.gt.0) then
        read(INTRAJ(m+1:m+k-1),*) Jintra(i)
       elseif(INTRAJ(m+1:m+9).ne.bl) then
        read(INTRAJ(m+1:m+9),*) J; Jintra(i) = 2*J
       else
       if(Jshell(i).eq.0) Jintra(i) = Jintra(i-1)
       if(Jintra(i-1).eq.0) Jintra(i) = Jshell(i)
       end if

      End do

      if(k.gt.0) then
       read(INTRAJ(m+1:m+k-1),*) Jintra(no)
      else
       k=INDEX(INTRAJ(m+1:m+9),'+')
       if(k.eq.0) k=INDEX(INTRAJ(m+1:m+9),'-')
       read(INTRAJ(m+1:m+k-1),*) J; Jintra(no) = 2*J
      end if

      Do i=1,no
       m = jn(i)+1
       if(Vshell(i).ne.0.or.iq(i).eq.m) Cycle
       if(iq(i).eq.1.or.iq(i).eq.m-1) then
         Vshell(i) = 1
       elseif(iq(i).eq.2.or.iq(i).eq.m-2) then
         if(Jshell(i).gt.0) Vshell(i) = 2
       elseif(iq(i).eq.3) then
         Vshell(i) = 3; if(Jshell(i).eq.jn(i)) Vshell(i)=1
       end if
      End do

      Jtotal = Jintra(no)

      End Subroutine Decode_cj


!======================================================================
      Subroutine Decode_confj
!======================================================================
!     decodes configuration from c-file format to INTEGER format
!----------------------------------------------------------------------
      Use conf_jj
      Implicit none
      Integer :: i,k,m

      no = INDEX(CONFIG,')',BACK=.TRUE.) / 9
      m=0
      Do i=1,no
       Call EL_nljk(CONFIG(m+1:m+5),nn(i),kn(i),ln(i),jn(i),in(i))
       read(CONFIG(m+7:m+8),'(i2)') iq(i)
       k = (2*ln(i)-jn(i))*(jn(i)+1)/2
       np_symc(i) = k*1000 + iq(i)
       m = m + 9
      End do

      End Subroutine Decode_confj

!======================================================================
      Subroutine Incode_confj
!======================================================================
!     incodes the configuration from INTEGER format to c-file format
!----------------------------------------------------------------------
      Use conf_jj
      Implicit none
      Integer :: i,m
      Character(5), External :: ELi

      m=0
      Do i=1,no
       if(iq(i).le.0) Cycle
       CONFIG(m+1:m+5) = ELi(nn(i),kn(i),in(i))
       write(CONFIG(m+6:m+9),'(a1,i2,a1)') '(',iq(i),')'
       m = m + 9
      End do

      ia = no*9

      End Subroutine Incode_confj


!======================================================================
      Subroutine Incode_confj1
!======================================================================
!     incodes the configuration (with index 1 in module conf_jj)
!     from INTEGER format to c-file format
!----------------------------------------------------------------------
      Use conf_jj
      Implicit none
      Integer :: i,m
      Character(5), External :: ELi

      m=0
      Do i=1,no1
       if(iq(i).le.0) Cycle
       CONFIG(m+1:m+5) = ELi(nn1(i),kn1(i),in1(i))
       write(CONFIG(m+6:m+9),'(a1,i2,a1)') '(',iq1(i),')'
       m = m + 9
      End do

      ia = no1*9

      End Subroutine Incode_confj1


!======================================================================
      Subroutine Decode_cjj(CONFIG,SHELLJ,INTRAJ,no,nn,kn,ln,jn,iq,in,&
                            Jshell,Vshell,Jintra)
!======================================================================
!     decode the configuration from c-file format to INTEGER format
!     Call: EL_nljk
!----------------------------------------------------------------------

      Implicit none
      Character(*), Intent(in) :: CONFIG,SHELLJ,INTRAJ
      Integer, Intent(out) :: no,nn(*),kn(*),ln(*),jn(*),iq(*),in(*),&
                              Jshell(*),Vshell(*),Jintra(*)
      Integer :: i,j,k,m

      m=INDEX(CONFIG,')',BACK=.TRUE.); no=m/9

      Vshell(1:no)=0; Jshell(1:no)=0; Jintra(1:no)=0
      m=-9
      Do i=1,no
       m = m + 9

       Call EL_nljk(CONFIG(m+1:m+5),nn(i),kn(i),ln(i),jn(i),in(i))
       read(CONFIG(m+7:m+8),'(i2)') iq(i)

       if(SHELLJ(m+5:m+5).eq.';') then
        read(SHELLJ(m+4:m+4),'(i1)') Vshell(i)
        read(SHELLJ(m+9:m+9),'(i1)') J; Jshell(i) = 2*J
       elseif(SHELLJ(m+8:m+8).eq.'/') then
        read(SHELLJ(m+1:m+7),'(i7)') Jshell(i)
       elseif(SHELLJ(m+9:m+9).ne.' ') then
        read(SHELLJ(m+1:m+9),'(i9)') J; Jshell(i) = 2*J
       end if

       k = INDEX(INTRAJ(m+1:m+9),'/')
       if(i.eq.1) then
        Jintra(i) = Jshell(i)
       elseif(i.eq.no) then
        Cycle
       elseif(k.gt.0) then
        read(INTRAJ(m+1:m+k-1),*) Jintra(i)
       elseif(LEN_TRIM(INTRAJ(m+1:m+9)).ne.0) then
        read(INTRAJ(m+1:m+9),*) J; Jintra(i) = 2*J
       else
       if(Jshell(i).eq.0) Jintra(i) = Jintra(i-1)
       if(Jintra(i-1).eq.0) Jintra(i) = Jshell(i)
       end if

      End do

      if(k.gt.0) then
       read(INTRAJ(m+1:m+k-1),*) Jintra(no)
      else
       k=INDEX(INTRAJ(m+1:m+9),'+')
       if(k.eq.0) k=INDEX(INTRAJ(m+1:m+9),'-')
       read(INTRAJ(m+1:m+k-1),*) J; Jintra(no) = 2*J
      end if

      Do i=1,no
       m = jn(i)+1
       if(Vshell(i).ne.0.or.iq(i).eq.m) Cycle
       if(iq(i).eq.1.or.iq(i).eq.m-1) then
         Vshell(i) = 1
       elseif(iq(i).eq.2.or.iq(i).eq.m-2) then
         if(Jshell(i).gt.0) Vshell(i) = 2
       elseif(iq(i).eq.3) then
         Vshell(i) = 3; if(Jshell(i).eq.jn(i)) Vshell(i)=1
       end if
      End do

      End Subroutine Decode_cjj


!======================================================================
      Subroutine Decode_configuration(conf)
!======================================================================
!     decodes configuration from the string
!----------------------------------------------------------------------
      Use conf_jj
      Implicit none
      Integer :: i,i1,i2,start
      Character(*) :: conf
      Character(5) :: EL, ELi

      no = 0
      start = 1
      Do
       i = index(conf(start:),')')
       if(i.eq.0) Exit
       no = no +1
       start = start+i
      End do

      if(no.eq.0) Return

      Call Clean_a(conf)
      start = 1
      Do i=1,no
       i1 = index(conf(start:),'(')+start-1
       i2 = index(conf(start:),')')+start-1
       EL = conf(start:i1-1)
       read(conf(i1+1:i2-1),*) iq(i)
       Call EL_NLJK(EL,nn(i),kn(i),ln(i),jn(i),in(i))
       start = i2+1
      End do

      End Subroutine Decode_configuration


!======================================================================
      Subroutine Decode_core(string)
!======================================================================
!     decodes core subshells  from the string
!----------------------------------------------------------------------
      Use conf_jj
      Implicit none
      Integer :: i,ii
      Character(*) :: string
      Character(5), external :: ELi

      Do i = 1,mcore
       read(string,*,IOSTAT=ii) e_core(1:i)
       if (ii /= 0) Exit
      End do
      ncore = i-1

      Do i = 1,ncore
       e_core(i) = adjustl(e_core(i))
       Call EL_NLJK(e_core(i),nn_core(i),k_core(i),l_core(i),j_core(i),ii)
       e_core(i) = ELi(nn_core(i),k_core(i),0)
      End do

      End Subroutine Decode_core


!======================================================================
      Subroutine Reduce_LS_jj(no,nn,kn,ln,jn,iq,in,iq_min,iq_max)
!======================================================================
! ... convert LS- to jj-configuration
!----------------------------------------------------------------------
      Implicit none
      Integer :: no,i,n,l,k,j,jj
      Integer :: nn(*),kn(*),ln(*),jn(*),iq(*),in(*),iq_min(*),iq_max(*)
      Integer :: n1(no),l1(no),q1(no),i1(no)

      n1(1:no) = nn(1:no)
      l1(1:no) = ln(1:no)
      q1(1:no) = iq(1:no)
      i1(1:no) = in(1:no)

      n=0
      Do i=1,no; l=l1(i)

       if(l.eq.0) then
        j=l+l+1; k=(l+l-j)*(j+1)/2
        n = n + 1
        nn(n) = n1(i); kn(n)=k; ln(n)=l; jn(n)=j; iq(n)=0; in(n)=i1(i)
        iq_min(n)=q1(i); iq_max(n)=q1(i)
        Cycle
       end if

       j=l+l-1; k=(l+l-j)*(j+1)/2; jj=l+l+1
       n = n + 1
       nn(n) = n1(i); kn(n)=k; ln(n)=l; jn(n)=j; iq(n)=0; in(n)=i1(i)
       iq_min(n)=max(0,q1(i)-jj-1); iq_max(n)=min(q1(i),j+1)

       j=l+l+1; k=(l+l-j)*(j+1)/2; jj=l+l-1
       n = n + 1
       nn(n) = n1(i); kn(n)=k; ln(n)=l; jn(n)=j; iq(n)=q1(i); in(n)=i1(i)
       iq_min(n)=max(0,q1(i)-jj-1); iq_max(n)=min(q1(i),j+1)

      End do
      no = n

      End Subroutine Reduce_LS_jj


!======================================================================
      Subroutine Reduce_jj_LS(no,nn,kn,ln,jn,iq,in)
!======================================================================
! ... convert jj- to LS-configuration
!----------------------------------------------------------------------
      Implicit none
      Integer :: no,i,n
      Integer :: nn(*),kn(*),ln(*),jn(*),iq(*),in(*)

      n=1
      Do i=2,no
       if(nn(n).eq.nn(i).and.ln(n).eq.ln(i)) then
        iq(n)=iq(n)+iq(i)
       else
        n=n+1
        nn(n)=nn(i); ln(n)=ln(i); iq(n)=iq(i); in(n)=in(i)
       end if
      End do
      no = n

      End Subroutine Reduce_jj_LS
