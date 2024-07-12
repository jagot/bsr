set(BSR_RECOUP4_COMMON_SOURCES
  MOD_bsr_recoup.f90
  add_block.f90
  matrix_updates.f90
  read_arg.f90)

add_executable(bsr_recoup4
  ${BSR_RECOUP4_COMMON_SOURCES}
  bsr_recoup.f90
  record_matrix.f90
  sub1.f90)
target_link_libraries(bsr_recoup4 zconfls zcom shells)

if(${BUILD_MPI_PROGRAMS})
  add_executable(bsr_recoup4_mpi
    ${BSR_RECOUP4_COMMON_SOURCES}
    bsr_recoup_mpi.f90
    record_matrix_mpi.f90
    sub1_mpi.f90
    br_arg.f90
    br_channel.f90
    br_channels_ion.f90
    br_target.f90
    br_target_ion.f90
    collect_ACF.f90)
  target_link_libraries(bsr_recoup4_mpi zconfls zcom shells MPI::MPI_Fortran)
endif()
