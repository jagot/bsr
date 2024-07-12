if(${BUILD_MPI_PROGRAMS})
  set_source_files_properties(crees.f PROPERTIES COMPILE_OPTIONS "-std=legacy")
  add_executable(bsr_phot3_mpi
    mod_phot.f90
    ak_coef.f90
    bsr_phot.f90
    phot_sec.f90
    read_arg.f90
    read_Ddat.f90
    read_Hdat.f90
    read_Wdat.f90
    sort_photo.f90
    sym_mat1.f90
    zaface.f90
    zopen.f90
    read_bound.f90
    crees.f)
  target_link_libraries(bsr_phot3_mpi zcom ${BLAS_LAPACK_LIBRARIES} MPI::MPI_Fortran)
endif()
