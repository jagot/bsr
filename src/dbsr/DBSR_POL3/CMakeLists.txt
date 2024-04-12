add_executable(dbsr_pol3
  mod_pol.f90              dbsr_pol.f90             add_nortb.f90
  read_data.f90            r_arg.f90                r_dbsrmat.f90
  r_dipmat.f90             solv_mat.f90)
target_link_libraries(dbsr_pol3 dbs zconfjj zcom ${BLAS_LAPACK_LIBRARIES})
