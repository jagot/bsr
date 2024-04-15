add_executable(dbsr_ci3
  mod_dbsr_ci.f90    read_arg.f90       setup_matrix.f90
        add_integral.f90   dbsr_ci.f90        diag.f90
        gen_matrix.f90     int_corr.f90       inf_dbsr_ci.f90
        read_conf_jj.f90   read_data.f90      read_shift.f90
        rk_evaluate.f90    sk_evaluate.f90)
target_link_libraries(dbsr_ci3 dbs zconfjj zcom ${BLAS_LAPACK_LIBRARIES})
