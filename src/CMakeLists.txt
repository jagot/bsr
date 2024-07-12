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

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set (Fortran_RELEASE_FLAGS  "-funroll-all-loops -fno-automatic -fcray-pointer -fno-f2c -O2 -march=native -fallow-argument-mismatch")
  set (Fortran_DEBUG_FLAGS    "-fno-f2c -O0 -g -fcheck=all -fbacktrace")
  set (Fortran_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")
else (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (Fortran_RELEASE_FLAGS "-O2")
  set (Fortran_DEBUG_FLAGS   "-O0 -g")
  set (Fortran_COVERAGE_FLAGS "")
endif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

if (CMAKE_C_COMPILER_ID MATCHES "GNU")
  set (C_RELEASE_FLAGS  "-O3 -march=native -Wall -Wextra")
  set (C_DEBUG_FLAGS    "-O0 -g -Wall -Wextra")
  set (C_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")
elseif (CMAKE_C_COMPILER_ID MATCHES "Clang|AppleClang")
  set (C_RELEASE_FLAGS "-O3 -march=native -Wall -Wextra")
  set (C_DEBUG_FLAGS   "-O0 -g -Wall -Wextra")
  set (C_COVERAGE_FLAGS "-fprofile-instr-generate -fcoverage-mapping")
else (CMAKE_C_COMPILER_ID MATCHES "GNU")
  message ("CMAKE_C_COMPILER full path: " ${CMAKE_C_COMPILER})
  message ("C compiler: " ${CMAKE_C_COMPILER_ID})
  message ("No optimized C compiler flags are known, we just try -O2...")
  set (C_RELEASE_FLAGS "-O2")
  set (C_DEBUG_FLAGS   "-O0 -g -DDEBUG")
  set (C_COVERAGE_FLAGS "")
endif (CMAKE_C_COMPILER_ID MATCHES "GNU")

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set (CXX_RELEASE_FLAGS  "-O3 --std=c++11 -march=native")
  set (CXX_DEBUG_FLAGS    "-O0 -g")
  set (CXX_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang")
  set (CXX_RELEASE_FLAGS "-O3 --std=c++11 -march=native")
  set (CXX_DEBUG_FLAGS   "-O0 -g")
  set (CXX_COVERAGE_FLAGS "-fprofile-instr-generate -fcoverage-mapping")
else (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  message ("CMAKE_CXX_COMPILER full path: " ${CMAKE_CXX_COMPILER})
  message ("C compiler: " ${CMAKE_CXX_COMPILER_ID})
  message ("No optimized C compiler flags are known, we just try -O2...")
  set (CXX_RELEASE_FLAGS "-O2")
  set (CXX_DEBUG_FLAGS   "-O0 -g")
  set (CXX_COVERAGE_FLAGS "")
endif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")

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

find_package(BLAS)
find_package(LAPACK)
find_package(MPI)
if(${MPI_Fortran_FOUND} AND ${BUILD_MPI_PROGRAMS})
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

set(FC_MALLOC LINUX)
set(BLAS_LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

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
