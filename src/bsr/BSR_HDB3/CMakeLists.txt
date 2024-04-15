add_executable(bsr_hdb3
        mod_blacs.f90            mod_hd.f90               add_exp.f90
        add_mvc.f90              bsr_hdb.f90              b_out.f90
        diag_hd.f90              h_out.f90                pri_mainpar.f90
        read_diag.f90            rsol_out.f90             read_arg.f90
        sub1_hd.f90              transform_mat.f90        transform_ovl.f90
        w_out.f90                find_channel_label.f90   find_channel_nodes.f90)
message("${SCALAPACK_LIBRARIES}")
target_link_libraries(bsr_hdb3 bs zconfls zcom scalapack ${BLAS_LAPACK_LIBRARIES} MPI::MPI_Fortran)
