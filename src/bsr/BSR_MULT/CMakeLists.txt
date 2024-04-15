if(BSR_VERSION STREQUAL "4")
  add_subdirectory(mult4)
  add_subdirectory(mult4_mpi)
else()
  add_subdirectory(mult3)
  add_subdirectory(mult3_tab)
  add_subdirectory(simple2)
endif()

