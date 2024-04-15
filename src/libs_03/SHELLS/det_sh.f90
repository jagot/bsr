!====================================================================
      Subroutine  DET_sh  (l,q,id,ML,MS,Idet)
!====================================================================
!
!     gives the determinant id and ML,MS for the subshell lq
!     Input:  lq - subshell
!             id - the determinant pointer
!     Output: Idet - description of the determinant as
!                    list of orbital positions in the common list
!             ML,MS - total azimutal quantum numbers (2*m+1)
!     Calls:  ML_id, MS_id
!
!---------------------------------------------------------------------

      Use det_exp; Use det_f3; Use det_f4; Use det_f5

      Implicit none
      Integer,Intent(in ) :: l,q,id
      Integer,Intent(out) :: ML,MS, Idet(q)
      Integer :: i,j,k,nd,nq
      Integer, External :: ML_id, MS_id

!----------------------------------------------------------------------
      if(id.le.0.or.l.lt.0) Call DET_shw (l,q,id)
!----------------------------------------------------------------------
! ... two-electron case for l>5, the only case we should compute:

      if(q.eq.2.and.l.gt.5) then
       nq=4*l+2; nd=nq*(nq-1)/2; if(id.gt.nd) Call DET_shw (l,q,id)
       k=0
       Do i=1,nq-1
        Do j=i+1,nq
         k=k+1; if(k.eq.id) Exit
        End do; if(k.eq.id) Exit
       End do
       Idet(1)=i; Idet(2)=j
       ML=ML_id(i)+ML_id(j)-1; MS=MS_id(i)+MS_id(j)-1
       Return

!----------------------------------------------------------------------

      elseif(q.eq.1) then                    ! l^1 - subshell:
       if(id.gt.4*l+2) Call DET_shw (l,q,id)
       ML=ML_id(id); MS=MS_id(id); Idet(1)=id; Return

      elseif(q.eq.4*l+2) then                ! full subshell:
       if(id.gt.1) Call DET_shw (l,q,id)
       ML=1; MS=1; Do i=1,q; Idet(i)=i; End do; Return

      elseif(q.eq.4*l+1) then                ! one-electron vacancy:
       if(id.gt.4*l+2) Call DET_shw (l,q,id)
       ML=-ML_id(id)+2; MS=-MS_id(id)+2
       k=0
       Do i=1,4*l+2
        if(i.ne.id) then; k=k+1; Idet(k)=i;  end if
       End do
       Return

      else

!----------------------------------------------------------------
      SELECT CASE(l*100+q)

      CASE(102)                                  !  p2 - subshell
       if(id.gt.kd_p2) Call DET_shw (l,q,id)
       ML=ML_p2(id); MS=MS_p2(id); Idet(1:q)= Idet_p2(1:q,id)
      CASE(103)                                  !  p3 - subshell
       if(id.gt.kd_p3) Call DET_shw (l,q,id)
       ML=ML_p3(id); MS=MS_p3(id); Idet(1:q)= Idet_p3(1:q,id)
      CASE(104)                                  !  p4 - subshell
       if(id.gt.kd_p4) Call DET_shw (l,q,id)
       ML=ML_p4(id); MS=MS_p4(id); Idet(1:q)= Idet_p4(1:q,id)
      CASE(202)                                  !  d2 - subshell
       if(id.gt.kd_d2) Call DET_shw (l,q,id)
       ML=ML_d2(id); MS=MS_d2(id); Idet(1:q)= Idet_d2(1:q,id)
      CASE(203)                                  !  d3 - subshell
       if(id.gt.kd_d3) Call DET_shw (l,q,id)
       ML=ML_d3(id); MS=MS_d3(id); Idet(1:q)= Idet_d3(1:q,id)
      CASE(204)                                  !  d4 - subshell
       if(id.gt.kd_d4) Call DET_shw (l,q,id)
       ML=ML_d4(id); MS=MS_d4(id); Idet(1:q)= Idet_d4(1:q,id)
      CASE(205)                                  !  d5 - subshell
       if(id.gt.kd_d5) Call DET_shw (l,q,id)
       ML=ML_d5(id); MS=MS_d5(id); Idet(1:q)= Idet_d5(1:q,id)
      CASE(206)                                  !  d6 - subshell
       if(id.gt.kd_d6) Call DET_shw (l,q,id)
       ML=ML_d6(id); MS=MS_d6(id); Idet(1:q)= Idet_d6(1:q,id)
      CASE(207)                                  !  d7 - subshell
       if(id.gt.kd_d7) Call DET_shw (l,q,id)
       ML=ML_d7(id); MS=MS_d7(id); Idet(1:q)= Idet_d7(1:q,id)
      CASE(208)                                  !  d8 - subshell
       if(id.gt.kd_d8) Call DET_shw (l,q,id)
       ML=ML_d8(id); MS=MS_d8(id); Idet(1:q)= Idet_d8(1:q,id)
      CASE(302)                                  !  f2 - subshell
       if(id.gt.kd_f2) Call DET_shw (l,q,id)
       ML=ML_f2(id); MS=MS_f2(id); Idet(1:q)= Idet_f2(1:q,id)
      CASE(312)                                  !  f12 - subshell
       if(id.gt.kd_f12) Call DET_shw (l,q,id)
       ML=ML_f12(id); MS=MS_f12(id); Idet(1:q)= Idet_f12(1:q,id)
      CASE(402)                                  !  g2 - subshell
       if(id.gt.kd_g2) Call DET_shw (l,q,id)
       ML=ML_g2(id); MS=MS_g2(id); Idet(1:q)= Idet_g2(1:q,id)
      CASE(502)                                  !  h2 - subshell
       if(id.gt.kd_h2) Call DET_shw (l,q,id)
       ML=ML_h2(id); MS=MS_h2(id); Idet(1:q)= Idet_h2(1:q,id)

      CASE(303)                                  !  f3 - subshell
       if(id.gt.kd_f3) Call DET_shw (l,q,id)
       ML=ML_f3(id); MS=MS_f3(id); Idet(1:q)= Idet_f3(1:q,id)
      CASE(304)                                  !  f4 - subshell
       if(id.gt.kd_f4) Call DET_shw (l,q,id)
       ML=ML_f4(id); MS=MS_f4(id); Idet(1:q)= Idet_f4(1:q,id)
      CASE(305)                                  !  f3 - subshell
       if(id.gt.kd_f5) Call DET_shw (l,q,id)
       ML=ML_f5(id); MS=MS_f5(id); Idet(1:q)= Idet_f5(1:q,id)

      Case default
        Call DET_shw (l,q,id)
        Stop ' DET_sh: lq - out of possible values'
      END SELECT
      end if

      End  Subroutine  DET_sh


!----------------------------------------------------------------------
      Subroutine DET_shw (l,q,id)
!----------------------------------------------------------------------

      Integer, Intent(in) :: l,q,id

      write(6,'(a,a,3i5)') 'DET_sh: id is out the range: ', &
                           'l,q,id = ',l,q,id
      Stop  'DET_sh: it(id) is out the range:'

      End Subroutine DET_shw
