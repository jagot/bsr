add_executable(dbsr_hf
  dbsr_lib_zcom.f90            dbsr_lib_zconfjj.f90      dbsr_lib_dbs.f90
  hf_MOD_dbsr_hf.f90           hf_MOD_df_orbitals.f90    hf_MOD_energy.f90
  dbsr_hf.f90
  hf_boundary.f90              hf_check_tails.f90        hf_def_conf_jj.f90
  hf_def_conf_LS.f90           hf_energy_BR.f90          hf_energy_SE.f90
  hf_energy_VP.f90             hf_get_case.f90           hf_get_estimates.f90
  hf_get_spline_param.f90      hf_matrix.f90             hf_orthogonality.f90
  hf_plot_bsw.f90              hf_read_pqbs.f90          hf_read_GRASP.f90
  hf_rk.f90                    hf_rotate.f90             hf_sk.f90
  hf_solve_HF.f90              hf_state_energies.f90     hf_summry.f90
  write_nl.f90)
target_link_libraries(dbsr_hf ${BLAS_LAPACK_LIBRARIES})
