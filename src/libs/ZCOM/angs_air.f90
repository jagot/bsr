!======================================================================
      Function ANGS_AIR (angs)
!======================================================================
!     correction to the wavelength in air - taken from CFF MCHF package
!----------------------------------------------------------------------
      Use zconst

      Implicit none
      Real(dp), Intent(in) :: angs
      Real(dp) :: sigma, ANGS_AIR

      ANGS_AIR=angs
      if (ANGS .GT. 2000e0_dp) then
        SIGMA = (1.e8_dp/ANGS)**2
        ANGS_AIR = ANGS/(one  + 8342.13e-8_dp                &
                              + 2406030/(130.e+8_dp - SIGMA) &
                              + 15997  /(38.9e+8_dp - SIGMA) )
      end if

      End Function ANGS_AIR