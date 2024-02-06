!======================================================================
      Subroutine Conv_au (Z,AWT,au_cm,au_eV,ipri)
!======================================================================
!     "au_cm" is transformation from au to cm^-1.
!     It depends on the atomic weight AWT according to:
!        au_cm = 2*Ry_cm*AWT/(AWT+RATIO) cm-1/au,
!     where RATIO - the electron mass in au
!           Ry_cm - Rydberg constatnt in cm-1
!           Ry_eV - Rydberg constatnt in eV
!     These constants are taken from  http://physics.nist.gov (2001)
!----------------------------------------------------------------------
      Implicit none
      Real(8), intent(in) :: Z,AWT
      Real(8), intent(out) :: au_cm, au_ev
      Integer, intent(in) :: ipri
      Real(8) :: A
      Real(8), parameter :: RATIO = 5.485799110d-4
      Real(8), parameter :: Ry_cm = 109737.31568549d0
      Real(8), parameter :: Ry_ev = 13.6056917d0

! ... determine suitable atomic weight if any

      if(AWT.ne.0.d0) then
        A = AWT
      else
        if ( Z .eq. 1.) then
           A = 1.
         else if ( Z .gt. 10.) then
           A = 2*Z+1 + (Z-11)/2
         else if ( MOD(INT(Z),2) .eq. 0 .or. Z .eq. 7. ) then
           A = 2*Z
         else
           A = 2*Z+1
        end if
      end if

      au_cm = 2*Ry_cm * A / (A + RATIO)
      au_eV = 2*Ry_ev * A / (A + RATIO)

      if(ipri.gt.0) then
       write(ipri,'(/a,f10.2)') 'Atomic weight = ',A
       write(ipri,'(a,f20.10)') &
            'a.u. --> cm-1 , au_cm = ', au_cm
       write(ipri,'(a,f20.10/)') &
            'a.u. --> eV   , au_eV = ', au_eV
      end if

      End Subroutine Conv_au
