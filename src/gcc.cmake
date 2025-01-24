set (Fortran_RELEASE_FLAGS  "-O3 -march=native -fallow-argument-mismatch -g1 -floop-block -fexternal-blas -fblas-matmul-limit=50 -cpp -fopenmp")
set (Fortran_DEBUG_FLAGS    "-O0 -g -fcheck=all -fbacktrace -cpp -fopenmp")
set (Fortran_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage -cpp -fopenmp")

set (C_RELEASE_FLAGS  "-O3 -march=native -Wall -Wextra")
set (C_DEBUG_FLAGS    "-O0 -g -Wall -Wextra")
set (C_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")

set (CXX_RELEASE_FLAGS  "-O3 --std=c++11 -march=native")
set (CXX_DEBUG_FLAGS    "-O0 -g")
set (CXX_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")

include("${CMAKE_CURRENT_LIST_DIR}/speedups.cmake")
