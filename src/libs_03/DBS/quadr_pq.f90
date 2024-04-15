!======================================================================
      Real(8) Function quadr_pq (io,jo,m)
!======================================================================
!     Evaluates   <p_io | r^m | p_jo>     with respect to r
!     between orbitals io and jo
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: io,jo,m
      Integer :: i,j, iv,ith,jth
      Real(8) :: pi,pj,qi,qj

      quadr_pq = 0.d0

      Do iv=1,nv
       gx(:) = gr(iv,:)**m * grw(iv,:)
       Do ith=1,ksp;  i=iv+ith-1; pi=pq(i,1,io)
       gw(:) = gx(:)*pbsp(iv,:,ith)
       Do jth=1,ksp;  j=iv+jth-1; pj=pq(j,1,jo)
         quadr_pq = quadr_pq + SUM(pbsp(iv,:,jth)*gw(:))*pi*pj
       End do
       End do
      End do

      Do iv=1,nv
       gx(:) = gr(iv,:)**m * grw(iv,:)
       Do ith=1,ksq;  i=iv+ith-1; qi=pq(i,2,io)
       gw(:) = gx(:)*qbsp(iv,:,ith)
       Do jth=1,ksq;  j=iv+jth-1; qj=pq(j,2,jo)
         quadr_pq = quadr_pq + SUM(qbsp(iv,:,jth)*gw(:))*qi*qj
       End do
       End do
      End do

      End Function quadr_pq


!======================================================================
      Real(8) Function quadr_qp (io,jo,m,ip,jp)
!======================================================================
!     Evaluates   <p_io | r^m | p_jo>     with respect to r
!     where p  - two-component Dirac functions in module DF_orbitals
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_orbitals_pq

      Implicit none
      Integer, intent(in) :: io,jo,m,ip,jp
      Real(8), external :: quadrm

      Select Case (ip*10+jp)

      Case(11)
       quadr_qp = quadrm(nv,ks,ksp,ksp,m,pbsp,pbsp,pq(1,1,io),pq(1,1,jo))
      Case(12)
       quadr_qp = quadrm(nv,ks,ksp,ksq,m,pbsp,qbsp,pq(1,1,io),pq(1,2,jo))
      Case(21)
       quadr_qp = quadrm(nv,ks,ksq,ksp,m,qbsp,pbsp,pq(1,2,io),pq(1,1,jo))
      Case(22)
       quadr_qp = quadrm(nv,ks,ksq,ksq,m,qbsp,qbsp,pq(1,2,io),pq(1,2,jo))

      Case default
       Stop 'quadr_qp: inknown combination ip,jp'

      End Select

      End Function quadr_qp


!======================================================================
      Real(8) Function quadr_00 (io,jo)
!======================================================================
!     Evaluates overlaps  <p_io | p_jo>    between orbitals io and jo
!----------------------------------------------------------------------
      Use DBS_orbitals_pq

      quadr_00 = SUM(pq(:,:,io)*bpq(:,:,jo))

      End Function quadr_00


