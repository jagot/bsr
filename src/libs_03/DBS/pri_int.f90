!=======================================================================
      Subroutine pri_int (nu,icase,k,i1,i2,i3,i4,S,SS)
!=======================================================================
!     print the integral; idexes are given for orbitals
!-----------------------------------------------------------------------
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: nu,icase,k,i1,i2,i3,i4
      Real(8), intent(in) :: S,SS
      Character AINT(0:3)/'O','L','R','S'/

      if(icase.eq.1) then
       write(nu,'(a,a,a,a,a,a,2D20.8)') &
        AINT(icase),' (',EBS(i1),',',EBS(i2),') = ',S,SS
      else
       write(nu,'(a,i2,a,a,a,a,a,a,a,a,a,D20.8,D15.5)') &
        AINT(icase),k,' (',EBS(i1),',',EBS(i2),';', &
        EBS(i3),',',EBS(i4),') = ',S,SS
      end if

      End Subroutine pri_int


!=======================================================================
      Subroutine prj_int (nu,icase,k,i1,i2,i3,i4,S,SS,ic,jc,atype)
!=======================================================================
!     debuging printing the integral values
!-----------------------------------------------------------------------
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: nu,icase,k,i1,i2,i3,i4,ic,jc,atype
      Real(8), intent(in) :: S,SS
      Character AINT(0:3)/'O','L','R','S'/
      Character(1) :: P1,P2,P3,P4

      if(icase.eq.1) then
       write(nu,'(a,a,a,a,a,a,2D20.10)') &
        AINT(icase),' (',EBS(i1),',',EBS(i2),') = ',S,SS

      elseif(icase.eq.2) then
        Select case(atype)
         Case(0);  P1='P'; P2='P'; P3='P'; P4='P'
         Case(1);  P1='Q'; P2='Q'; P3='Q'; P4='Q'
         Case(2);  P1='P'; P2='Q'; P3='P'; P4='Q'
         Case(3);  P1='Q'; P2='P'; P3='Q'; P4='P'
        End Select

       write(nu,'(T14,a,i2,13a,E15.7)') &
        AINT(icase),k,' (',P1,ebs(i1),',',P2,ebs(i2),';', &
        P3,EBS(i3),',',P4,EBS(i4),') = ',S

      elseif(icase.eq.3) then
        Select case(atype)
         Case(0);  P1='P'; P2='P'; P3='Q'; P4='Q'
         Case(1);  P1='P'; P2='Q'; P3='Q'; P4='P'
        End Select

       write(nu,'(a,i2,13a,2F10.5,2i6)') &
        AINT(icase),k,' (',P1,ebs(i1),',',P2,ebs(i2),';', &
        P3,EBS(i3),',',P4,EBS(i4),') = ',S,SS,ic,jc
      end if

      End Subroutine prj_int

