!======================================================================
      Subroutine convol(ns,ks,a,d,icase,sym_i,sym_j)
!======================================================================
!     convolutes the rkb(i,j,i',j') array of spline integrals
!     with density matrix d(:,:) 
!
!     results in array a(:,:)
!
!     icase =  1  - convolution other 2 and 4 variables, RK(.a;.b)
!              2  - convolution other 1 and 3 variables, RK(a.;b.)
!              3  - convolution other 2 and 3 variables, RK(.a;b.)
!              4  - convolution other 1 and 4 variables, RK(a.;.b)
!
!     sym_i  ->  symmetry in respect of i,i' variables ('s','l','n')
!     sym_j  ->  symmetry in respect of j,j' variables ('s','l','n')
!
!     combination of sym_i and sym_j leads to different represantation
!     for a and d:   a(ns,ks),  a(ns,2*ks+1),  a(ns,ns)
!----------------------------------------------------------------------
      Use DBS_integrals, only: rkb
      Use DBS_debug

      Implicit none
      Integer, intent(in) :: ns,ks,icase
      Character, intent(in) :: sym_i,sym_j
      Real(8), intent(in ) :: d(ns,*)
      Real(8), intent(out) :: a(ns,*)
      ! local variables
      Integer :: i,j, ip,jp, imin,imax, jmin,jmax, ii,jj
      Real(8) :: c,t1,t2

      t1 = RRTC()

      if(icase.le.2) a(1:ns,1:ks)=0.d0
      if(icase.gt.2) a(1:ns,1:ns)=0.d0

      Select case(icase)
!----------------------------------------------------------------------
      Case(1)                                          !  I( . a ; . b)

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

        do ip=1,ks
        do i=1,ns-ip+1   
           a(i,ip) = SUM(d(1:ns,1:ks)*rkb(i,1:ns,ip,1:ks))
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=0.d0
        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do ip=1,ks
        do i=ks+1-ip,ns;  c=0.d0
        do jp=1,ks
        do j=ks+1-jp,ns;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do ip=1,ks
        do i=ks+1-ip,ns;  c=0.d0
        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=0.d0
        do jp=1,ks
        do j=ks+1-jp,ns;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      end if
!----------------------------------------------------------------------
      Case(2)                                         !  I( a . ; b . )

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

        do ip=1,ks
        do i=1,ns-ip+1   
           a(i,ip) = SUM(d(1:ns,1:ks)*rkb(1:ns,i,1:ks,ip))
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=0.d0
        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns;  c=0.d0
        do ip=1,ks
        do i=ks+1-ip,ns;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=0.d0
        do ip=1,ks
        do i=ks+1-ip,ns;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns;  c=0.d0
        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      end if 

!----------------------------------------------------------------------
      Case(3);  a(1:ns,1:ns) = 0.d0                    ! I( . a ; b . )

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

       do jp=1,ks
       do j =1,ns-jp+1; jj=j+jp-1
       do ip=1,ks
       do i =1,ns-ip+1; ii=i+ip-1
          c = rkb(i,j,ip,jp);      a( i,jj) = a( i,jj) + c*d(ii, j)
          if(ip.gt.1)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(jp.gt.1)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(ip.gt.1.and.jp.gt.1)  a(ii, j) = a(ii, j) + c*d( i,jj)
       end do;  end do
       end do;  end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks 
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks 
          c = rkb(i,j,ip,jp);       a( i,jj) = a( i,jj) + c*d(ii, j)
          if(ip.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(jp.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(ip.ne.ks.and.jp.ne.ks) a(ii, j) = a(ii, j) + c*d( i,jj)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

       do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax; jj=j+jp-ks
       do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
       do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);      a( i,jj) = a( i,jj) + c*d(ii, j)
       end do; end do
       end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks 
        do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);      a( i,jj) = a( i,jj) + c*d(ii, j)
          if(jp.ne.ks)             a( i, j) = a( i, j) + c*d(ii,jj)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks 
          c = rkb(i,j,ip,jp);       a( i,jj) = a( i,jj) + c*d(ii, j)
          if(ip.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
        end do; end do
        end do; end do

      end if
!----------------------------------------------------------------------
      Case(4);  a(1:ns,1:ns) = 0.d0                    ! I( a . ; . b )

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

       do jp=1,ks
       do j =1,ns-jp+1; jj=j+jp-1
       do ip=1,ks
       do i =1,ns-ip+1; ii=i+ip-1
          c = rkb(i,j,ip,jp);      a(ii, j) = a(ii, j) + c*d( i,jj)
          if(ip.gt.1)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(jp.gt.1)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(ip.gt.1.and.jp.gt.1)  a( i,jj) = a( i,jj) + c*d(ii, j)
       end do;  end do
       end do;  end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

       do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax; jj=j+jp-ks
       do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
       do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);      a(ii, j) = a(ii, j) + c*d( i,jj)
       end do; end do
       end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks 
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks 
          c = rkb(i,j,ip,jp);       a(ii, j) = a(ii, j) + c*d( i,jj)
          if(ip.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(jp.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(ip.ne.ks.and.jp.ne.ks) a( i,jj) = a( i,jj) + c*d(ii, j)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks 
        do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a(ii, j) = a(ii, j) + c*d( i,jj)
          if(jp.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks 
          c = rkb(i,j,ip,jp);       a(ii, j) = a(ii, j) + c*d( i,jj)
          if(ip.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
        end do; end do
        end do; end do

      end if
!----------------------------------------------------------------------
      Case default

       Stop 'convol: unknown case'

      End Select

      t2 = RRTC()
      ic_convol = ic_convol + 1
      if(icase.le.2) time_convol=time_convol + (t2-t1)

      End Subroutine convol
