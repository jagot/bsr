add_subdirectory(BSR_BREIT3)
add_subdirectory(BSR_BREIT4)

if(${MPI_Fortran_FOUND})
  add_subdirectory(BSR_BREIT3_MPI)
  add_subdirectory(BSR_BREIT4_MPI)
endif()
