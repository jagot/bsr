if(${SCALAPACK_FOUND})
   add_executable(dbsr_hd3_mpi
     mod_blacs.f90            mod_hd.f90               add_exp.f90
     b_out.f90                dbsr_hdb.f90             diag_mat.f90
     h_out.f90                check_dbsr_mat.f90       read_overlaps.f90
     read_matrix.f90          rsol_out.f90             read_arg.f90
     sub1_hd.f90              w_out.f90                channel_label.f90)
   target_link_libraries(dbsr_hd3_mpi dbs zconfjj zcom scalapack ${BLAS_LAPACK_LIBRARIES} MPI::MPI_Fortran)
endif()
