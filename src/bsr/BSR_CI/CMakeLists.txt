add_executable(bsr_ci
     MOD_bsr_ci.f90  MOD_c_data.f90    MOD_l_core.f90  MOD_z_core.f90
     bsr_ci.f90      read_arg.f90      read_conf.f90   read_bsw.f90
     matrix_LS.f90   diag.f90          matrix_rel.f90  matrix_J.f90
     inform.f90      read_int_bnk.f90  I_data.f90      g_factor.f90
     read_data.f90)
target_link_libraries(bsr_ci bs zconfls zcom ${BLAS_LAPACK_LIBRARIES})
