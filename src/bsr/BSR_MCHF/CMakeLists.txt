add_executable(bsr_mchf
     MOD_bsr_mchf.f90   bsr_mchf.f90
     get_case.f90       read_conf.f90     get_spline_param.f90
     read_ang_coef.f90  def_orbitals.f90  boundary.f90
     get_estimates.f90  check_tails.f90   diag.f90
     update_hs.f90      mchf_matrix.f90   ab.f90
     rotate.f90         scf_mchf.f90      solve_eiv.f90
     solve_direct.f90   solve_nr.f90      orthogonality.f90
     write_bsw.f90      summry.f90        ecore_hf.f90
     hlc_hf.f90)
target_link_libraries(bsr_mchf bs zconfls zcom ${BLAS_LAPACK_LIBRARIES})
