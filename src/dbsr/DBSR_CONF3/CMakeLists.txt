set(DBSR_CONF_COMMON_SOURCES
  mod_dbsr_conf.f90     gen_conf.f90       make_coupling.f90
  check_perturber.f90   def_phys_targ.f90  read_arg.f90      sub_jj.f90
  record_orth.f90       def_orth_cond.f90  def_del.f90       check_cfg.f90
  dbsr_conf_inf.f90)

add_executable(dbsr_conf3 ${DBSR_CONF_COMMON_SOURCES} dbsr_conf.f90)
target_link_libraries(dbsr_conf3 zconfjj zcom)

if(${MPI_Fortran_FOUND})
   add_executable(dbsr_conf3_mpi ${DBSR_CONF_COMMON_SOURCES}
    dbsr_conf_mpi.f90   br_arg.f90  br_conf_jj.f90  br_phys_orb_jj.f90  br_target_jj.f90)
   target_link_libraries(dbsr_conf3_mpi zconfjj zcom bsr_mpi)
endif()
