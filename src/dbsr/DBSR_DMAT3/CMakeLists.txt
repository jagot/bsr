add_executable(dbsr_dmat3
  mod_dbsr_dmat.f90        mod_cmdata.f90           check_det.f90
  check_mult_bnk.f90       det_fact.f90             dbsr_dmat.f90
  dv_out.f90               d_data.f90               d_matr.f90
  d_out.f90                dd_out.f90
  gen_matrix.f90           gen_zf.f90
  get_Vdip.f90             idef_type.f90            gen_dbs.f90
  read_arg.f90             read_data.f90            updates.f90)
target_link_libraries(dbsr_dmat3 dbs zconfjj zcom)
