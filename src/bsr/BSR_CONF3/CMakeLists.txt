set(BSR_CONF_COMMON_SOURCES
  mod_bsr_conf.f90
  check_cfg.f90
  check_pertuber.f90
  def_orth_cond.f90
  def_phys_targ.f90
  gen_conf.f90
  make_coupling.f90
  record_orth.f90
  sub_jj.f90
  sub_jk.f90
  sub_ls.f90
  read_arg.f90
  bsr_conf_inf.f90)

add_executable(bsr_conf3
  ${BSR_CONF_COMMON_SOURCES}
  bsr_conf.f90)
target_link_libraries(bsr_conf3 zconfls zcom shells)

if(${BUILD_MPI_PROGRAMS})
  add_executable(bsr_conf3_mpi
    ${BSR_CONF_COMMON_SOURCES}
    bsr_conf_mpi.f90
    br_arg.f90)
  target_link_libraries(bsr_conf3_mpi bsr_mpi zconfls zcom shells MPI::MPI_Fortran)
endif()
