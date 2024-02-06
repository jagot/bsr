!======================================================================
      Subroutine convol(ns,ks,a,d,icase,sym_i,sym_j)
!======================================================================
!
!    convolutes the rkb(i,j,i',j') array of spline integrals
!    with density matrix d(:,:)
!
!    results in array a
!
!    icase =  1  - convolution other 2 and 4 variables, RK(.a;.b)
!             2  - convolution other 1 and 3 variables, RK(a.;b.)
!             3  - convolution other 2 and 3 variables, RK(.a;b.)
!             4  - convolution other 1 and 4 variables, RK(a.;.b)
!
!    sym_i  ->  symmetry in respect of i,i' variables ('s', 'n')
!    sym_j  ->  symmetry in respect of j,j' variables ('s', 'n')
!
!----------------------------------------------------------------------
      Use spline_integrals, only: rkb

      Implicit none
      Integer, intent(in) :: ns,ks,icase
      Character, intent(in) :: sym_i,sym_j
      Real(8), intent(in ) :: d(ns,*)
      Real(8), intent(out) :: a(ns,*)

      Integer :: i,j, ip,jp, imin,imax, jmin,jmax, ii,jj
      Real(8) :: c


      if(icase.gt.2) a(1:ns,1:ns) = 0.d0

      if(sym_i.eq.'s'.and.sym_j.eq.'s'.and.icase.eq.1) then

        do ip = 1,ks
          do i = 1,ns-ip+1
            c = 0.d0
            do jp = 1,ks
              do j = 1,ns-jp+1
                c = c + d(j,jp)*rkb(i,j,ip,jp)
              end do
            end do
            a(i,ip) = c
          end do
        end do

      elseif(sym_i.eq.'s'.and.sym_j.eq.'s'.and.icase.eq.2) then

        do ip = 1,ks
          do i = 1,ns-ip+1
            c = 0.d0
            do jp = 1,ks
              do j = 1,ns-jp+1
                c = c + d(j,jp)*rkb(j,i,jp,ip)
              end do
            end do
            a(i,ip) = c
          end do
        end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n'.and.icase.eq.1) then

        do ip = 1,ks+ks-1
          imin=max( 1, 1 + ks-ip)
          imax=min(ns,ns + ks-ip)
          do i = imin,imax
            c = 0.d0
            do jp = 1,ks+ks-1
              jmin=max( 1, 1 + ks-jp)
              jmax=min(ns,ns + ks-jp)
              do j = jmin,jmax
                  c = c + d(j,jp)*rkb(i,j,ip,jp)
              end do
            end do
            a(i,ip) = c
          end do
        end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n'.and.icase.eq.2) then

        do ip = 1,ks+ks-1
          imin=max( 1, 1 + ks-ip)
          imax=min(ns,ns + ks-ip)
          do i = imin,imax
            c = 0.d0
            do jp = 1,ks+ks-1
              jmin=max( 1, 1 + ks-jp)
              jmax=min(ns,ns + ks-jp)
              do j = jmin,jmax
                  c = c + d(j,jp)*rkb(j,i,jp,ip)
              end do
            end do
            a(i,ip) = c
          end do
        end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'s'.and.icase.eq.1) then

        do ip = 1,ks+ks-1
          imin=max( 1, 1 + ks-ip)
          imax=min(ns,ns + ks-ip)
          do i = imin,imax
            c = 0.d0
            do jp = 1,ks
              do j = 1,ns-jp+1
                  c = c + d(j,jp)*rkb(i,j,ip,jp)
              end do
            end do
            a(i,ip) = c
          end do
        end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'s'.and.icase.eq.2) then

        do jp = 1,ks
          do j = 1,ns-jp+1
            c = 0.d0
            do ip = 1,ks+ks-1
              imin=max( 1, 1 + ks-ip)
              imax=min(ns,ns + ks-ip)
              do i = imin,imax
                  c = c + d(i,ip)*rkb(i,j,ip,jp)
              end do
            end do
            a(j,jp) = c
          end do
        end do


      elseif(sym_i.eq.'s'.and.sym_j.eq.'n'.and.icase.eq.2) then

        do jp = 1,ks+ks-1
          jmin=max( 1, 1 + ks-jp)
          jmax=min(ns,ns + ks-jp)
          do j = jmin,jmax
            c = 0.d0
            do ip = 1,ks
              do i = 1,ns-ip+1
                  c = c + d(i,ip)*rkb(i,j,ip,jp)
              end do
            end do
            a(j,jp) = c
          end do
        end do

      elseif(sym_i.eq.'s'.and.sym_j.eq.'n'.and.icase.eq.1) then

        do ip = 1,ks
          do i = 1,ns-ip+1
            c = 0.d0
            do jp = 1,ks+ks-1
              jmin=max( 1, 1 + ks-jp)
              jmax=min(ns,ns + ks-jp)
              do j = jmin,jmax
                  c = c + d(j,jp)*rkb(i,j,ip,jp)
              end do
            end do
            a(i,ip) = c
          end do
        end do


    elseif(sym_i.eq.'s'.and.sym_j.eq.'s'.and.icase.eq.3) then

      do jp = 1,ks
        do j = 1,ns-jp+1
          jj = j+jp-1
          do ip = 1,ks
            do i = 1,ns-ip+1
              ii = i+ip-1

              c = rkb(i,j,ip,jp)

                                a( i,jj) = a( i,jj) + c*d(ii, j)
       if(ip.gt.1)              a(ii,jj) = a(ii,jj) + c*d( i, j)
       if(jp.gt.1)              a( i, j) = a( i, j) + c*d(ii,jj)
       if(ip.gt.1.and.jp.gt.1)  a(ii, j) = a(ii, j) + c*d( i,jj)

            end do
          end do
        end do
      end do

    elseif(sym_i.eq.'s'.and.sym_j.eq.'s'.and.icase.eq.4) then

      do jp = 1,ks
        do j = 1,ns-jp+1
          jj = j+jp-1
          do ip = 1,ks
            do i = 1,ns-ip+1
              ii = i+ip-1

              c = rkb(i,j,ip,jp)

                                a(ii, j) = a(ii, j) + c*d( i,jj)
       if(ip.gt.1)              a( i, j) = a( i, j) + c*d(ii,jj)
       if(jp.gt.1)              a(ii,jj) = a(ii,jj) + c*d( i, j)
       if(ip.gt.1.and.jp.gt.1)  a( i,jj) = a( i,jj) + c*d(ii, j)

            end do
          end do
        end do
      end do

    elseif(sym_i.eq.'n'.and.sym_j.eq.'n'.and.icase.eq.3) then

      do jp = 1,ks+ks-1
        jmin=max( 1, 1 + ks-jp)
        jmax=min(ns,ns + ks-jp)
        do j = jmin,jmax
          jj = j+jp-ks
          do ip = 1,ks+ks-1
            imin=max( 1, 1 + ks-ip)
            imax=min(ns,ns + ks-ip)
            do i = imin,imax
            ii = i+ip-ks

            a( i,jj) = a( i,jj) + rkb(i,j,ip,jp)*d(ii, j)

            end do
          end do
        end do
      end do

    elseif(sym_i.eq.'n'.and.sym_j.eq.'n'.and.icase.eq.4) then

      do jp = 1,ks+ks-1
        jmin=max( 1, 1 + ks-jp)
        jmax=min(ns,ns + ks-jp)
        do j = jmin,jmax
          jj = j+jp-ks
          do ip = 1,ks+ks-1
            imin=max( 1, 1 + ks-ip)
            imax=min(ns,ns + ks-ip)
            do i = imin,imax
            ii = i+ip-ks

            a(ii, j) = a(ii, j) + rkb(i,j,ip,jp)*d( i,jj)

            end do
          end do
        end do
      end do

    elseif(sym_i.eq.'n'.and.sym_j.eq.'s'.and.icase.eq.3) then

      do jp = 1,ks
        do j = 1,ns-jp+1
          jj = j+jp-1
          do ip = 1,ks+ks-1
            imin=max( 1, 1 + ks-ip)
            imax=min(ns,ns + ks-ip)
            do i = imin,imax
              ii = i+ip-ks

                     a( i,jj) = a( i,jj) + rkb(i,j,ip,jp)*d(ii, j)
         if(jp.gt.1) a( i, j) = a( i, j) + rkb(i,j,ip,jp)*d(ii,jj)

            end do
          end do
        end do
      end do


      elseif(sym_i.eq.'n'.and.sym_j.eq.'s'.and.icase.eq.4) then

      do jp = 1,ks
        do j = 1,ns-jp+1
          jj = j+jp-1
          do ip = 1,ks+ks-1
            imin=max( 1, 1 + ks-ip)
            imax=min(ns,ns + ks-ip)
            do i = imin,imax
              ii = i+ip-ks

                     a(ii, j) = a(ii, j) + rkb(i,j,ip,jp)*d( i,jj)
         if(jp.gt.1) a(ii,jj) = a(ii,jj) + rkb(i,j,ip,jp)*d( i, j)

            end do
          end do
        end do
      end do


    elseif(sym_i.eq.'s'.and.sym_j.eq.'n'.and.icase.eq.3) then

      do jp = 1,ks+ks-1
        jmin=max( 1, 1 + ks-jp)
        jmax=min(ns,ns + ks-jp)
        do j = jmin,jmax
          jj = j+jp-ks
          do ip = 1,ks
            do i = 1,ns-ip+1
              ii = i+ip-1

                     a( i,jj) = a( i,jj) + rkb(i,j,ip,jp)*d(ii, j)   ! ????
       if(ip.gt.1)   a(ii,jj) = a(ii,jj) + rkb(i,j,ip,jp)*d( i, j)

            end do
          end do
        end do
      end do

     elseif(sym_i.eq.'s'.and.sym_j.eq.'n'.and.icase.eq.4) then

      do jp = 1,ks+ks-1
        jmin=max( 1, 1 + ks-jp)
        jmax=min(ns,ns + ks-jp)
        do j = jmin,jmax
          jj = j+jp-ks
          do ip = 1,ks
            do i = 1,ns-ip+1
              ii = i+ip-1

                     a(ii, j) = a(ii, j) + rkb(i,j,ip,jp)*d( i,jj)   ! ????
       if(ip.gt.1)   a( i, j) = a( i, j) + rkb(i,j,ip,jp)*d(ii,jj)

            end do
          end do
        end do
      end do


      else

        Stop ' CONVOL:  unknown case '

      end if

      End Subroutine convol
