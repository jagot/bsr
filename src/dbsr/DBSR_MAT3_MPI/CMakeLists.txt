set(DBSR_MAT3A_SOURCE_COMMON
  mod_dbsr_mat.f90
  L_data.f90
  O_data.f90
  R_data.f90
  S_data.f90
  add_integral.f90
  check_mat.f90
  dbs_orth.f90
  diag.f90
  gen_matrix.f90
  idef_type.f90
  read_arg.f90
  read_data.f90
  state_res.f90           
  target_subs.f90
  updates.f90
  zero_cond.f90)

add_executable(dbsr_mat3a
  ${DBSR_MAT3A_SOURCE_COMMON}
  br_fake.f90
  dbsr_mat.f90
  record_matrix.f90
  sub1.f90)
target_link_libraries(dbsr_mat3a dbs zconfjj zcom ${BLAS_LAPACK_LIBRARIES})

if(${BUILD_MPI_PROGRAMS})
  add_executable(dbsr_mat3a_mpi
    ${DBSR_MAT3A_SOURCE_COMMON}
    dbsr_mat_mpi.f90
    MOD_c_data.f90
    record_matrix_mpi.f90
    sub1_mpi.f90
    check_det_new.f90
    det_fact_new.f90
    br_DBS_orbitals_pq.f90
    br_arg.f90
    br_buffer.f90
    br_channel_jj.f90
    br_conf_jj.f90
    br_dets.f90
    br_orb_jj.f90
    br_grid_jj.f90
    br_phys_orb_jj.f90
    br_radial_overlaps.f90
    br_target_jj.f90
    br_barrier.f90)
  target_link_libraries(dbsr_mat3a_mpi dbs zconfjj zcom
    ${BLAS_LAPACK_LIBRARIES} MPI::MPI_Fortran)
endif()
