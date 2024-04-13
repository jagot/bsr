set(BSR_MAT4_COMMON_SOURCES
  MOD_bsr_mat.f90          MOD_c_data.f90           l_core.f90               z_core.f90
  bs_orth.f90              check_det.f90            f_values.f90
  get_v.f90                gen_matrix.f90           get_orth_chan.f90
  idef_type.f90            jsym_int.f90
  I_data.f90               L_data.f90               O_data.f90
  read_arg.f90             read_bsw.f90             matrix_updates.f90
  read_conf.f90            read_data.f90            rw_mask.f90
  Z_data.f90               check_idef.f90           bsr_mat_inf.f90)

add_executable(bsr_mat4
  ${BSR_MAT4_COMMON_SOURCES}
  bsr_mat.f90
  read_matrix.f90
  record_matrix.f90
  state_res.f90
  sub1.f90
  check_mat.f90)
target_link_libraries(bsr_mat4 bs zconfls zcom ${BLAS_LAPACK_LIBRARIES})

add_executable(bsr_mat4_mpi
  ${BSR_MAT4_COMMON_SOURCES}
  bsr_mat_mpi.f90
  br_arg.f90
  check_mat_mpi.f90
  collect_mpi.f90
  read_matrix_mpi.f90
  record_matrix_mpi.f90
  state_res_mpi.f90
  stop_mpi.f90
  sub1_mpi.f90)
target_link_libraries(bsr_mat4_mpi bs zconfls zcom bsr_mpi MPI::MPI_Fortran ${BLAS_LAPACK_LIBRARIES})
