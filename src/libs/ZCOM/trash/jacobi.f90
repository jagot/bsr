!======================================================================
      SUBROUTINE Jacobi(a,n,np,d,v,nrot)
!======================================================================
!     Computes all eigenvalues and eigenvectors of a real symmetric
!     matrix a, which is of size n by n, stored in a physical np by np
!     array.
!
!     On output, elements of a above the diagonal are destroyed.
!     d returns the eigenvalues of a  in its first n elements.
!     v is a matrix with the same logical and physical dimensions as a,
!     whose columns contain, on output, the normalized eigenvectors of a.
!     nrot returns the number of Jacobi rotations that were required.
!
!     after Numerical Recipes, f77
!----------------------------------------------------------------------
      Integer :: n,np,nrot
      Real(8) :: a(np,np),d(np),v(np,np)

      Integer :: i,j,ip,iq
      Real(8) :: c,g,h,s,sm,t,tau,theta,tresh,b(n),z(n)

! ... Initialize to the identity matrix.

      do ip=1,n
       do iq=1,n
        v(ip,iq)=0.d0
       end do
       v(ip,ip)=1.d0
      end do

! ... Initialize b and d to the diagonal of a.
! ... Vector z will accumulate terms of the form t*a_pq
! ... as in equation (11.1.14).

      do ip=1,n
       b(ip)=a(ip,ip)
       d(ip)=b(ip)
       z(ip)=0.d0
      end do
      nrot=0

      do i=1,50    !  sweeps

! ... Sum off-diagonal elements.

      sm=0.
      do ip=1,n-1
       do iq=ip+1,n
        sm=sm+abs(a(ip,iq))
       end do
      end do

! ... The normal return, which relies on quadratic convergence
! ... to machine underflow.

      if(sm.eq.0.)return

      if(i.lt.4)then
       tresh=0.2d0*sm/n**2       ! on the first three sweeps
      else
       tresh=0.d0                ! thereafter
      endif

      do ip=1,n-1
      do iq=ip+1,n
       g=100.d0*abs(a(ip,iq))

! ...  After four sweeps, skip the rotation
! ...  if the off-diagonal element is small.

      if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and. &
         (abs(d(iq))+g.eq.abs(d(iq)))) then
       a(ip,iq)=0.d0
      elseif(abs(a(ip,iq)).gt.tresh) then
       h=d(iq)-d(ip)
       if(abs(h)+g.eq.abs(h))then
        t=a(ip,iq)/h                ! t=1/(2*theta)
       else
        theta=0.5d0*h/a(ip,iq)      ! Equation (11.1.10)
        t=1.d0/(abs(theta)+sqrt(1.+theta**2))
        if(theta.lt.0.d0) t=-t
       end if

       c=1.d0/sqrt(1.d0+t**2)
       s=t*c
       tau=s/(1.d0+c)
       h=t*a(ip,iq)
       z(ip)=z(ip)-h
       z(iq)=z(iq)+h
       d(ip)=d(ip)-h
       d(iq)=d(iq)+h
       a(ip,iq)=0.d0

       do j=1,ip-1             ! Case of rotations 1 <= j < p
        g=a(j,ip)
        h=a(j,iq)
        a(j,ip)=g-s*(h+g*tau)
        a(j,iq)=h+s*(g-h*tau)
       end do

       do j=ip+1,iq-1          ! Case of rotations p < j < q
        g=a(ip,j)
        h=a(j,iq)
        a(ip,j)=g-s*(h+g*tau)
        a(j,iq)=h+s*(g-h*tau)
       end do

       do j=iq+1,n             ! Case of rotations q < j <= n
        g=a(ip,j)
        h=a(iq,j)
        a(ip,j)=g-s*(h+g*tau)
        a(iq,j)=h+s*(g-h*tau)
       end do

       do j=1,n
        g=v(j,ip)
        h=v(j,iq)
        v(j,ip)=g-s*(h+g*tau)
        v(j,iq)=h+s*(g-h*tau)
       end do

       nrot=nrot+1

      endif

      end do ! over iq
      end do ! over ip

! ... Update d with the sum of t*a_pq, and reinitialize z.

      do ip=1,n
       b(ip)=b(ip)+z(ip)
       d(ip)=b(ip)
       z(ip)=0.d0
      end do

      end do ! over sweeps i

      Stop 'jacobi: too many iterations'

      END SUBROUTINE jacobi

