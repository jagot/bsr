!======================================================================
      Subroutine density(ns,ks,d,p1,p2,type)
!======================================================================
!     d - density of two vectors p1,p2
!     type = 's' - symmetric upper-band column mode, or 'u'
!            'l' - symmetric lower-band column mode
!            'n' - non-symmetric band column mode
!            'x' - full matrix case, default
!----------------------------------------------------------------------
      Use DBS_debug

      Implicit none
      Character :: type
      Integer :: ns,ks, i,j,imin,imax
      Real(8) :: p1(ns),p2(ns)
      Real(8) :: d(ns,*)
      Real(8) :: t1,t2

      Call CPU_time(t1)

      if(type.eq.'s'.or.type.eq.'u') then

        d(1:ns,1:ks) = 0.d0
        do i=1,ns; d(i,1)=p1(i)*p2(i); end do
        do j=2,ks; do i=1,ns-j+1
          d(i,j)=p1(i)*p2(i+j-1) + p1(i+j-1)*p2(i)
        end do; end do

      elseif(type.eq.'l') then

        d(1:ns,1:ks) = 0.d0
        do i=1,ns; d(i,ks)= p1(i)*p2(i); end do
        do j=1,ks-1; do i=ks+1-j,ns
          d(i,j)=p1(i)*p2(i+j-ks) + p1(i+j-ks)*p2(i)
        end do; end do

      elseif(type.eq.'n') then

        d(1:ns,1:ks+ks-1) = 0.d0
        do j=1,ks+ks-1
          imin=max0( 1, 1+ks-j)
          imax=min0(ns,ns+ks-j)
          do i=imin,imax; d(i,j)=p1(i)*p2(i+j-ks); end do
        end do

      else

        d(1:ns,1:ns)=0.d0
        do i=1,ns; do j=1,ns
          d(i,j)=p1(i)*p2(j)
        end do; end do

      end if

      Call CPU_time(t2)
      ic_density=ic_density+1
      time_density = time_density + (t2-t1)

      End Subroutine density



