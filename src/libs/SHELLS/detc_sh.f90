!======================================================================
      Real(8)  Function  DETC_sh (l,q,it,id)
!======================================================================
!
!     gives the expansion coefficient for the given subshell 'lq',
!     term 'it' and determinant 'id'
!     Input:  lq - subshell
!             it - the term number
!             id - the determinant number
!     Calls:  CLEBSH, ML_id, MS_id
!
!----------------------------------------------------------------------

      Use det_exp; Use det_f3; Use det_f4; Use det_f5

      Implicit none
      Integer, Intent(in) :: l,q,it,id
      Integer :: i,j,k,ll,m,ml,ml1,ml2,ms,ms1,ms2,IL,IS,nd,nq
      Integer, External :: ML_id, MS_id
      Real(8) :: C, D1 = 1.d0
      Real(8), External :: Clebsh

!----------------------------------------------------------------------
      if(it.le.0.or.id.le.0) Call DETC_shw (l,q,it,id)
!----------------------------------------------------------------------

! ... two-electron case for l > 5, the only case we should compute:

      if(q.eq.2.and.l.gt.5) then
       if(it.gt.2*l+1) Call DETC_shw (l,q,it,id)
       IL=2*it-1; IS=2*mod(it-1,2)+1    ! 1S, 3P, 1D, 3F, ...
       nq=4*l+2; nd=nq*(nq-1)/2
       if(id.gt.nd) Call DETC_shw (l,q,it,id)
       k=0
       Do i=1,nq-1;
        Do j=i+1,nq
         k=k+1; if(k.eq.id) Exit
        End do; if(k.eq.id) Exit
       End do
       ml1=ML_id(i);  ml2=ML_id(j)
       ms1=MS_id(i);  ms2=MS_id(j)
       ML=ml1+ml2-1;  MS=ms1+ms2-1
       ll=l+l+1
       DETC_sh=CLEBSH(ll,ml1,ll,ml2,IL,ML)* &
               CLEBSH( 2,ms1, 2,ms2,IS,MS)*sqrt(2.d0)
        Return

!---------------------------------------------------------------------

      elseif(q.eq.1) then                    !  one-electron shell
       if(it.ne.1.or.id.gt.4*l+2) Call DETC_shw (l,q,it,id)
       DETC_sh = D1; Return

      elseif(q.eq.4*l+2) then                !  full subshell
       if(it.ne.1.or.id.ne.1) Call DETC_shw (l,q,it,id)
       DETC_sh = -D1; Return

      elseif(q.eq.4*l+1) then                !  one-electron vacancy
       if(it.ne.1.or.id.gt.4*l+2) Call DETC_shw (l,q,it,id)
       m=-ML_id(id)+2; m=iabs(m-1)/2
       DETC_sh=D1; if(mod(m+l,2).eq.1) DETC_sh=-D1; Return

      else
!----------------------------------------------------------------
      SELECT CASE(l*100+q)

      CASE(102)                                  !  p2 - subshell
       if(it.gt.nt_p2.or.id.gt.kd_p2) Call DETC_shw(l,q,it,id)
       i=Ip2(it,id); j=Np2(it)
      CASE(103)                                  !  p3 - subshell
       if(it.gt.nt_p3.or.id.gt.kd_p3) Call DETC_shw(l,q,it,id)
       i=Ip3(it,id); j=Np3(it)
      CASE(104)                                  !  p4 - subshell
       if(it.gt.nt_p4.or.id.gt.kd_p4) Call DETC_shw(l,q,it,id)
       i=Ip4(it,id); j=Np4(it)
      CASE(202)                                  !  d2 - subshell
       if(it.gt.nt_d2.or.id.gt.kd_d2) Call DETC_shw(l,q,it,id)
       i=Id2(it,id); j=Nd2(it)
      CASE(203)                                  !  d3 - subshell
       if(it.gt.nt_d3.or.id.gt.kd_d3) Call DETC_shw(l,q,it,id)
       i=Id3(it,id); j=Nd3(it)
      CASE(204)                                  !  d4 - subshell
       if(it.gt.nt_d4.or.id.gt.kd_d4) Call DETC_shw(l,q,it,id)
       i=Id4(it,id); j=Nd4(it)
      CASE(205)                                  !  d5 - subshell
       if(it.gt.nt_d5.or.id.gt.kd_d5) Call DETC_shw(l,q,it,id)
       i=Id5(it,id); j=Nd5(it)
      CASE(206)                                  !  d6 - subshell
       if(it.gt.nt_d6.or.id.gt.kd_d6) Call DETC_shw(l,q,it,id)
       i=Id6(it,id); j=Nd6(it)
      CASE(207)                                  !  d7 - subshell
       if(it.gt.nt_d7.or.id.gt.kd_d7) Call DETC_shw(l,q,it,id)
       i=Id7(it,id); j=Nd7(it)
      CASE(208)                                  !  d8 - subshell
       if(it.gt.nt_d8.or.id.gt.kd_d8) Call DETC_shw(l,q,it,id)
       i=Id8(it,id); j=Nd8(it)
      CASE(302)                                  !  f2 - subshell
       if(it.gt.nt_f2.or.id.gt.kd_f2) Call DETC_shw(l,q,it,id)
       i=If2(it,id); j=Nf2(it)
      CASE(312)                                  !  f12 - subshell
       if(it.gt.nt_f12.or.id.gt.kd_f12) Call DETC_shw(l,q,it,id)
       i=If12(it,id); j=Nf12(it)
      CASE(402)                                  !  g2 - subshell
       if(it.gt.nt_g2.or.id.gt.kd_g2) Call DETC_shw(l,q,it,id)
       i=Ig2(it,id); j=Ng2(it)
      CASE(502)                                  !  h2 - subshell
       if(it.gt.nt_h2.or.id.gt.kd_h2) Call DETC_shw(l,q,it,id)
       i=Ih2(it,id); j=Nh2(it)

      CASE(303)                                  !  f3 - subshell
       if(it.gt.nt_f3.or.id.gt.kd_f3) Call DETC_shw(l,q,it,id)
       i=ID_f3(id,it); j=JD_f3(id,it)
      CASE(304)                                  !  f4- subshell
       if(it.gt.nt_f4.or.id.gt.kd_f4) Call DETC_shw(l,q,it,id)
       i=ID_f4(id,it); j=JD_f4(id,it)
      CASE(305)                                  !  f5- subshell
       if(it.gt.nt_f5.or.id.gt.kd_f5) Call DETC_shw(l,q,it,id)
       i=ID_f5(id,it); j=JD_f5(id,it)

      CASE DEFAULT

       write(*,*) ' l,q = ',l,q
       Stop ' DETC_sh: lq - outside of possible values'

      END SELECT
      end if

      C=DBLE(i)/DBLE(j); C=sqrt(abs(C)); if(i.lt.0) C=-C; DETC_sh=C

      End  Function  DETC_sh

!----------------------------------------------------------------------
      Subroutine DETC_shw (l,q,it,id)
!----------------------------------------------------------------------

      Integer, Intent(in) :: l,q,it,id

	     write(6,'(a,a,4i5)') 'DETC_sh: it(id) is out the range: ', &
                           'l,q,it,id = ',l,q,it,id
      Stop  'DETC_sh: it(id) is out the range'

      End Subroutine DETC_shw
