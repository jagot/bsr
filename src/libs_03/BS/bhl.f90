!======================================================================
      Function BHL(i,j)
!======================================================================
!     <P_i| L |P_j>  with inclusion of rel.shift if rel = .true.
!----------------------------------------------------------------------
      Use spline_param
      Use spline_orbitals
      Use spline_hl
      Use spline_atomic

      Implicit none
      Integer, intent(in) :: i,j
      Real(8) :: BHL
      Real(8), external :: BVMV

      if (iabs(LBS(i)-LBS(j)) .NE. 0) Stop ' HL:  LI <> LJ'

      Call HLM (lbs(i))

      BHL = BVMV(ns,ks,hl,'a',pbs(1,i),pbs(1,j))

      End Function BHL

