!======================================================================
      Real(8) Function RRTC()
!======================================================================      
!     give the running time in seconds
!     (some old version to adjust for different compilers) 
!----------------------------------------------------------------------
      Implicit none
 
!     CHARACTER(LEN=8) :: D
!     CHARACTER(LEN=10) :: T
!     INTEGER :: id,ih,im,is,ims

      Real :: TM(2), ETIME !, DTIME
 
! ... Power station Fortran 4.0:

!     Call DATE_AND_TIME(date=D,time=T)
!     read(D,'(6x,i2)') id
!     read(T,'(3i2,1x,i3)') ih,im,is,ims
!     RRTC = id*86400 + ih*3600 + im*60 + is
!     RRTC = RRTC + ims/1000.d0

! ... Digital Fortran 6.0:

      RRTC = ETIME(TM); ! RRTC = TM(1)
 
      End Function RRTC
                                                                                        
