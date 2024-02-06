!======================================================================
      Subroutine Coul_wf(l,z,e,r,fc,gc,fcp,gcp,fcq,gcq, iout,ifail)
!======================================================================
!     Coulomb wavefunctions for different range of energies
!
!     l - orbital momentum
!     z - effective charge
!     e - energy
!     r - radial point
!     fc,gc - regular and irregular solutions
!     fcp, gcp - 1/k (fc)', 1/k (gc)'
!     fcq, gcq - (fcp)', (gcp)'
!
!======================================================================

      Implicit none

      Integer, intent(in) :: l,iout
      Integer, intent(out) :: ifail
      Real(8), intent(in) :: z,e,r
      Real(8), intent(out) :: fc,gc,fcp,gcp,fcq,gcq

      Real(8), parameter :: zero = 0.d0, one = 1.d0, two = 2.d0
      Integer :: ll
      Real(8) :: rho, eta, fl, k, kk, fac


       k   = dsqrt(dabs(e))
       eta = - z / k
       kk  = dsqrt(k)
       rho = r * k
       fl  = l
       ll  = l*(l+1)

       if ( e .ge. zero ) then

         Call Zcoulfg90 (rho,eta,fl, fc,gc,fcp,gcp, 0,ifail)

!         if (ifail .ne. 0) write (iout,'(2(a,i4),3(a,d16.8))') &
!                          ' Zcoulfg : ifail =',ifail,' l =',l,&
!                          ' rho =',rho,' eta =',eta


         ! ... normalize to delta on E ...

         fc=fc/kk; gc=gc/kk; fcp=fcp/kk; gcp=gcp/kk

         fac = (ll/rho + two*eta) / rho - one

         fcq = fac * fc; gcq = fac * gc

       else

        Stop 'coul_wf: E < 0.0'


       end if

       End Subroutine Coul_wf



