!====================================================================    
      MODULE spline_densities
!====================================================================
!     contains symmetric (lower-band column storage mode) 
!     and non-symmetric densities for two orbitals
!--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: idensity=0, ids,jds, idx,jdx, idq,jdq
      REAL(8), ALLOCATABLE :: ds(:,:), dx(:,:), dq(:,:)
      END MODULE spline_densities


!======================================================================
      Subroutine make_density(io,jo,sym)
!======================================================================
      Use spline_param
      Use spline_orbitals
      Use spline_densities

      Character :: sym
      Integer :: io,jo, i,j

      if(idensity.eq.0) then
       Allocate(ds(ns,ks), dx(ns,ns), dq(ns,ns))
       ds=0.d0; ids=0;jds=0; idx=0;jdx=0; idq=0;jdq=0
       idensity = 1
      end if

      if(io.lt.1.or.io.gt.nbf.or.jo.lt.1.or.jo.gt.nbf) &
         Stop 'make_density: orbital is out of range'

      if(sym.eq.'s') then
                                             
       if(ids.eq.io.and.jds.eq.jo) Return

       do i=1,ns;  ds(i,ks) = pbs(i,io)*pbs(i,jo); end do                                                                              !     o*
       do j = 1,ks-1                           
        do i = ks+1-j,ns
         ds(i,j) =  pbs(i,io)*pbs(i+j-ks,jo) + pbs(i+j-ks,io)*pbs(i,jo)
        end do
       end do

       ids = io; jds = jo

      elseif(sym.eq.'x') then

       if(idx.eq.io.and.jdx.eq.jo) Return

       do i = 1,ns
        do j = 1,ns
         dx(i,j) =  pbs(i,io)*pbs(j,jo)
        end do
       end do

       idx = io; jdx = jo

      elseif(sym.eq.'q') then

       if(idq.eq.io.and.jdq.eq.jo) Return

       do i = 1,ns
        do j = 1,ns
         dq(i,j) =  qbs(i,io)*qbs(j,jo)
        end do
       end do

       idq = io; jdq = jo

      end if

      End Subroutine make_density
