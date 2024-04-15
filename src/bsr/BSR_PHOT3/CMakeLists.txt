set_source_files_properties(crees.f PROPERTIES COMPILE_OPTIONS "-std=legacy")
add_executable(bsr_phot3
  mod_phot.f90
  ak_coef.f90
  bsr_phot.f90
  phot_sec.f90
  read_arg.f90
  read_bound.f90
  read_Ddat.f90
  read_Hdat.f90
  read_Wdat.f90
  sort_photo.f90
  zaface.f90
  crees.f)
target_link_libraries(bsr_phot3 zcom ${BLAS_LAPACK_LIBRARIES})
