if(DEBUG_SPEEDUPS)

# DEBUG_SPEEDUPS is the general compilation definition we use to debug
# speed-ups of various hot code paths, comparing the results with the
# previous slower implementation. Additionally, we have some
# specialized flags that enable debugging for code paths that are
# deemed correct, through some testing, and are hence disabled by
# default.
set(DEBUG_SPEEDUPS_FLAGS
"-DDEBUG_SPEEDUPS"
# "-DDEBUG_SPEEDUPS_CONVERT_PQ"
# "-DDEBUG_SPEEDUPS_RECORD_MATRIX"
# "-DDEBUG_SPEEDUPS_CONVOL_PAPB_SS"
# "-DDEBUG_SPEEDUPS_CONVOL_APBP_SS"
)
list(JOIN DEBUG_SPEEDUPS_FLAGS " " DEBUG_SPEEDUPS_FLAGS)

set (Fortran_RELEASE_FLAGS  "${Fortran_RELEASE_FLAGS} ${DEBUG_SPEEDUPS_FLAGS}")
set (Fortran_DEBUG_FLAGS    "${Fortran_DEBUG_FLAGS} ${DEBUG_SPEEDUPS_FLAGS}")
set (Fortran_COVERAGE_FLAGS "${Fortran_COVERAGE_FLAGS} ${DEBUG_SPEEDUPS_FLAGS}")
endif(DEBUG_SPEEDUPS)
