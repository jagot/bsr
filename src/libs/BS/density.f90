!======================================================================
      Subroutine density(ns,ks,d,p1,p2,type)
!======================================================================
!
!     d - density of two w.f. p1,p2
!     type = 's' - simmetrical direct case
!            'n' - nonsimmetrical direct
!            'x' - exchange case
!----------------------------------------------------------------------
      Implicit none

      INTEGER, INTENT(in) :: ns,ks
      REAL(8), DIMENSION(ns) :: p1,p2
      REAL(8), DIMENSION(ns,*) :: d
      CHARACTER(1) :: type
      INTEGER :: i,j,imin,imax


      if(type.eq.'s') then                    !     o***

        d(1:ns,1:ks) = 0.d0
                                              !     o***
        do i =1,ns                            !     o***
          d(i,1) =  p1(i)*p2(i)               !     o***
        end do                                !     o**
                                              !     o*
        do j = 2,ks                           !     o
          do i = 1,ns-j+1
           d(i,j) =  p1(i)*p2(i+j-1) + p1(i+j-1)*p2(i)
          end do
        end do

      elseif(type.eq.'n') then                !     o***
                                              !    *o***
        d(1:ns,1:ks+ks-1) = 0.d0


        do j = 1,ks+ks-1                      !   **o***
          imin=max0( 1, 1+ks-j)               !  ***o***
          imax=min0(ns,ns+ks-j)               !  ***o**
          do i = imin,imax                    !  ***o*
            d(i,j) = p1(i)*p2(i+j-ks)         !  ***o
          end do
        end do

      else

        d(1:ns,1:ns) = 0.d0

        do i = 1,ns
          do j = 1,ns
            d(i,j) =  p1(i)*p2(j)
          end do
        end do

      end if

      End Subroutine density
