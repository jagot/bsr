add_library(bsr_mpi
  br_bsorb.f90      br_channel.f90   br_conf_LS.f90  br_conf_LS1.f90
  br_core.f90       br_core1.f90     br_dets.f90     br_file.f90
  br_grid.f90       br_orb.f90       br_phys_orb.f90
  br_symc_LS.f90    br_symt_LS.f90   br_target.f90   br_radial_overlaps.f90)
target_link_libraries(bsr_mpi MPI::MPI_Fortran bs zconfls zcom)
set_target_properties(bsr_mpi PROPERTIES Fortran_MODULE_DIRECTORY ${GLOBAL_Fortran_MODULE_DIRECTORY})
