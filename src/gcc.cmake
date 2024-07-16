set (Fortran_RELEASE_FLAGS  "-funroll-all-loops -fno-automatic -fcray-pointer -fno-f2c -O2 -march=native -fallow-argument-mismatch")
set (Fortran_DEBUG_FLAGS    "-fno-f2c -O0 -g -fcheck=all -fbacktrace")
set (Fortran_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")

set (C_RELEASE_FLAGS  "-O3 -march=native -Wall -Wextra")
set (C_DEBUG_FLAGS    "-O0 -g -Wall -Wextra")
set (C_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")

set (CXX_RELEASE_FLAGS  "-O3 --std=c++11 -march=native")
set (CXX_DEBUG_FLAGS    "-O0 -g")
set (CXX_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")
