!====================================================================
      Subroutine  DET_sh_LS  (l,q,id,ML,MS,Idet)
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
      Integer, intent(in ) :: l,q,id
	     Integer, intent(out) :: ML, MS, Idet(q)
      Integer :: i,j,k,nd,nq
      Integer, External :: ML_id, MS_id

      if(l .lt.0) Call DET_LS_stop(l,q,id)
      if(q .le.0) Call DET_LS_stop(l,q,id)
      if(id.le.0) Call DET_LS_stop(l,q,id)

      if(q.eq.1) then                                  ! l^1 subshell
       if(id.gt.4*l+2) Call DET_LS_stop(l,q,id)
       ML=ML_id(id); MS=MS_id(id); Idet(1)=id

      elseif(q.eq.4*l+2) then                          ! full subshell
       if(id.gt.1) Call DET_LS_stop(l,q,id)
       ML=1; MS=1; Do i=1,q; Idet(i)=i; End do

      elseif(q.eq.4*l+1) then                          ! one-electron vacancy 
       if(id.gt.4*l+2) Call DET_LS_stop(l,q,id)
       ML=-ML_id(id)+2; MS=-MS_id(id)+2
       k=0
       Do i=1,4*l+2; if(i.eq.id) Cycle; k=k+1; Idet(k)=i; End do

      elseif(q.eq.2) then                              ! l^2 subshell

        Select case(l)
         case(1);           Call Det_p2(id,ML,MS,Idet)
         case(2);           Call Det_d2(id,ML,MS,Idet)
         case(3);           Call Det_f2(id,ML,MS,Idet)
         case(4);           Call Det_g2(id,ML,MS,Idet)
         case(5);           Call Det_h2(id,ML,MS,Idet)
         case default
         nq=4*l+2; nd=nq*(nq-1)/2; if(id.gt.nd) Call DET_LS_stop(l,q,id)
         k=0
         Do i=1,nq-1
         Do j=i+1,nq
          k=k+1; if(k.eq.id) Exit
         End do; if(k.eq.id) Exit
         End do
         Idet(1)=i; Idet(2)=j
         ML=ML_id(i)+ML_id(j)-1; MS=MS_id(i)+MS_id(j)-1
        End Select

      else             
        Select case(l)
         case(1);   
                       if(q.eq.3) Call Det_p3(id,ML,MS,Idet)
                       if(q.eq.4) Call Det_p4(id,ML,MS,Idet)
         case(2);      
                       if(q.eq.3) Call Det_d3(id,ML,MS,Idet)
                       if(q.eq.4) Call Det_d4(id,ML,MS,Idet)
                       if(q.eq.5) Call Det_d5(id,ML,MS,Idet)
                       if(q.eq.6) Call Det_d6(id,ML,MS,Idet)
                       if(q.eq.7) Call Det_d7(id,ML,MS,Idet)
                       if(q.eq.8) Call Det_d8(id,ML,MS,Idet)
         case(3);      
                       if(q.eq.3) Call Det_f3(id,ML,MS,Idet)
                       if(q.eq.4) Call Det_f4(id,ML,MS,Idet)
         case default
                       Call Det_LS_stop (l,q,id)
        End Select
       end if     

      End  Subroutine  DET_sh_LS


!======================================================================
      Subroutine DET_LS_stop (l,q,id)
!======================================================================       

      Integer, Intent(in) :: l,q,id

	     write(*,*) 'DET_sh_LS: input is out of range: l,q,id = ',l,q,id
      Stop  ' '

      End Subroutine DET_LS_stop