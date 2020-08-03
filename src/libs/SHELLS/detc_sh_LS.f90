!====================================================================
      Real(8) Function DETC_lq  (l,q,it,id)
!====================================================================
!
!     gives the determinant id and ML,MS for the subshell lq
!     Input:  lq - subshell
!             id - the determinant pointer
!     Output: Idet - description of the determinant as
!                    list of orbital positions in the common list
!             ML,MS - total azimutal quantum numbers (2*m+1)
!     Calls:  ML_id, MS_id,  DET_LS_stop, det_lq
!
!---------------------------------------------------------------------

      Implicit none 

      Integer, intent(in) :: l,q,it,id

      Integer :: i,j,k,nd,nq,m,ll, IL,IS, ML,ML1,ML2, MS,MS1,MS2
      Real(8) :: C, D1 = 1.d0
      Integer, External :: ML_id, MS_id
      Real(8), External :: Clebsh
      Real(8), External :: Detc_p2, Detc_p3, Detc_p4
      Real(8), External :: Detc_d2, Detc_d3, Detc_d4, Detc_d5, &
                           Detc_d6, Detc_d7, Detc_d8
      Real(8), External :: Detc_f2, Detc_f3, Detc_f4, Detc_f5, &
                           Detc_f6, Detc_f7, Detc_f8, Detc_f9, &
                           Detc_f10,Detc_f11, Detc_f12
      Real(8), External :: Detc_g2, Detc_h2

      if(l .lt.0) Call DETC_lq_stop(i,q,it,id)
      if(q .le.0) Call DETC_lq_stop(i,q,it,id)
      if(it.le.0) Call DETC_lq_stop(i,q,it,id)
      if(id.le.0) Call DETC_lq_stop(i,q,it,id)

      if(q.eq.1) then                                  ! l^1 subshell
       if(it.ne.1.or.id.gt.4*l+2) Call DETC_lq_stop(i,q,it,id)
       DETC_lq = D1

      elseif(q.eq.4*l+2) then                          ! full subshell
       if(it.ne.1.or.id.ne.1) Call DETC_lq_stop(i,q,it,id)
       DETC_lq = -D1

      elseif(q.eq.4*l+1) then                          ! one-electron vacancy 
       if(it.ne.1.or.id.gt.4*l+2) Call DETC_lq_stop(i,q,it,id)
       m=-ML_id(id)+2; m=iabs(m-1)/2
       DETC_lq = D1; if(mod(m+l,2).eq.1) DETC_lq = -D1

      elseif(q.eq.2) then                              ! l^2 subshell

        Select case(l)
         case(1);           DETC_lq = Detc_p2(it,id)
         case(2);           DETC_lq = Detc_d2(it,id)
         case(3);           DETC_lq = Detc_f2(it,id)
         case(4);           DETC_lq = Detc_g2(it,id)
         case(5);           DETC_lq = Detc_h2(it,id)
         case default
          if(it.gt.2*l+1) Call DETC_lq_stop (l,q,it,id)
          IL=2*it-1; IS=2*mod(it-1,2)+1    ! 1S, 3P, 1D, 3F, ...
          nq=4*l+2; nd=nq*(nq-1)/2 
          if(id.gt.nd) Call DETC_lq_stop (l,q,it,id)
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
          DETC_lq = CLEBSH(ll,ml1,ll,ml2,IL,ML)* &
                    CLEBSH( 2,ms1, 2,ms2,IS,MS)*sqrt(2.d0)
        End Select

      else             
        Select case(l)
         case(1);   
                       if(q.eq.3) DETC_lq = Detc_p3(it,id)
                       if(q.eq.4) DETC_lq = Detc_p4(it,id)
         case(2);      
                       if(q.eq.3) DETC_lq = Detc_d3(it,id)
                       if(q.eq.4) DETC_lq = Detc_d4(it,id)
                       if(q.eq.5) DETC_lq = Detc_d5(it,id)
                       if(q.eq.6) DETC_lq = Detc_d6(it,id)
                       if(q.eq.7) DETC_lq = Detc_d7(it,id)
                       if(q.eq.8) DETC_lq = Detc_d8(it,id)
         case(3);      
                       if(q.eq.3) DETC_lq = Detc_f3(it,id)
                       if(q.eq.4) DETC_lq = Detc_f4(it,id)
         case default
                       Call Detc_lq_stop (l,q,it,id)
        End Select
       end if     

      End  Function  DETC_lq


!======================================================================
      Subroutine DETC_lq_stop (l,q,it,id)
!======================================================================       

      Integer, Intent(in) :: l,q,it,id

	     write(*,*) 'DETC_lq_LS: input is out of range: l,q,it,id = ',l,q,it,id
      Stop  ' '

      End Subroutine DETC_lq_stop