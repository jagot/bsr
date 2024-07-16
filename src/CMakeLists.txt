cmake_minimum_required(VERSION 3.22)
project(bsr)

IF(NOT DEFINED CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build : None Debug Release Coverage."
       FORCE)
ENDIF(NOT DEFINED CMAKE_BUILD_TYPE)
message("Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(GLOBAL_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mods)
include_directories(${GLOBAL_Fortran_MODULE_DIRECTORY})

enable_language(Fortran)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  add_definitions("-DDEBUG")
endif()

set(CMAKE_Fortran_FLAGS_RELEASE "${Fortran_RELEASE_FLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG   "${Fortran_DEBUG_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${C_RELEASE_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG   "${C_DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CXX_RELEASE_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG   "${CXX_DEBUG_FLAGS}")

# https://stackoverflow.com/a/11437693
SET(CMAKE_Fortran_FLAGS_COVERAGE
    "${Fortran_DEBUG_FLAGS} ${Fortran_COVERAGE_FLAGS}"
    CACHE STRING "Flags used by the C++ compiler during coverage builds."
    FORCE )
SET(CMAKE_C_FLAGS_COVERAGE
    "${C_DEBUG_FLAGS} ${C_COVERAGE_FLAGS}"
    CACHE STRING "Flags used by the C compiler during coverage builds."
    FORCE )
SET(CMAKE_CXX_FLAGS_COVERAGE
    "${CXX_DEBUG_FLAGS} ${CXX_COVERAGE_FLAGS}"
    CACHE STRING "Flags used by the C++ compiler during coverage builds."
    FORCE )
SET(CMAKE_EXE_LINKER_FLAGS_COVERAGE
    ""
    CACHE STRING "Flags used for linking binaries during coverage builds."
    FORCE )
SET(CMAKE_SHARED_LINKER_FLAGS_COVERAGE
    ""
    CACHE STRING "Flags used by the shared libraries linker during coverage builds."
    FORCE )
MARK_AS_ADVANCED(
    CMAKE_Fortran_FLAGS_COVERAGE
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_EXE_LINKER_FLAGS_COVERAGE
    CMAKE_SHARED_LINKER_FLAGS_COVERAGE )

IF(NOT DEFINED BUILD_MPI_PROGRAMS)
  set (BUILD_MPI_PROGRAMS ON CACHE BOOL
       "Should MPI versions of the codes be builts? Default: ON"
       FORCE)
ENDIF(NOT DEFINED BUILD_MPI_PROGRAMS)

find_package(BLAS REQUIRED)
if(${BUILD_MPI_PROGRAMS})
  find_package(MPI)
  if(${MPI_Fortran_FOUND})
    find_package(SCALAPACK) # This defined an imported target
                            # "scalapack" which we can simply link to.
    message("MPI implementation found, enabling build of MPI executables")
    if(${SCALAPACK_FOUND})
      message("MPI and Scalapack found - enabling MPI Fortran")
      include_directories(${MPI_Fortran_INCLUDE_PATH})
    endif()
  else()
    message("Disabling build of MPI executables")
    set(BUILD_MPI_PROGRAMS OFF)
  endif()
else()
  message("Disabling build of MPI executables")
endif()

set(FC_MALLOC LINUX)
set(BLAS_LAPACK_LIBRARIES ${BLAS_LIBRARIES})

if (NOT BSR_VERSION)
  set(BSR_VERSION 3)
endif()

message("Will build version ${BSR_VERSION} of libraries and programs")

if(BSR_VERSION STREQUAL "4")
  add_subdirectory(libs)
else()
  add_subdirectory(libs_03)
endif()

add_subdirectory(bsr)
add_subdirectory(dbsr)
