add_executable(bsr_hf
     MOD_bsr_hf.f90       MOD_hf_orbitals.f90    MOD_hf_energy.f90     MOD_boef_list.f90
     bsr_hf.f90           coef_ee_1conf.f90      I_int.f90             rk.f90
     get_case.f90         get_spline_param.f90   get_estimates.f90     check_tails.f90
     quadr_hf.f90         solve_HF.f90           boundary.f90
     orthogonality.f90    bhl_hf.f90             rotate.f90            hf_matrix.f90
     update_hs.f90        write_cfg.f90          write_bsw.f90         write_nl.f90
     plot_bsw.f90         summry.f90             hf_nr.f90)
target_link_libraries(bsr_hf bs zconfls shells zcom ${BLAS_LAPACK_LIBRARIES})
