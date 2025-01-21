set (Fortran_RELEASE_FLAGS  "-funroll-all-loops -fno-automatic -fcray-pointer -fno-f2c -O3 -march=native -fallow-argument-mismatch -g1 -floop-block -fexternal-blas -fblas-matmul-limit=50 -cpp")
set (Fortran_DEBUG_FLAGS    "-fno-f2c -O0 -g -fcheck=all -fbacktrace -cpp")
set (Fortran_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage -cpp")

if(DEBUG_SPEEDUPS)
set (Fortran_RELEASE_FLAGS  "${Fortran_RELEASE_FLAGS} -DDEBUG_SPEEDUPS")
set (Fortran_DEBUG_FLAGS    "${Fortran_DEBUG_FLAGS} -DDEBUG_SPEEDUPS")
set (Fortran_COVERAGE_FLAGS "${Fortran_COVERAGE_FLAGS} -DDEBUG_SPEEDUPS")
endif(DEBUG_SPEEDUPS)

set (C_RELEASE_FLAGS  "-O3 -march=native -Wall -Wextra")
set (C_DEBUG_FLAGS    "-O0 -g -Wall -Wextra")
set (C_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")

set (CXX_RELEASE_FLAGS  "-O3 --std=c++11 -march=native")
set (CXX_DEBUG_FLAGS    "-O0 -g")
set (CXX_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")
