if(MPI_Fortran_FOUND)
message(WARNING "bsr_breit4_mpi does not build correctly")
# add_executable(bsr_breit4_mpi
#   MOD_bsr_breit.f90        bsr_breit_mpi.f90
#   add_it_oper.f90          add_res.f90              br_arg.f90
#   breit_me.f90             conf_calc_mpi.f90        conf_loop_mpi.f90
#   get_det_mpi.f90          get_res_mpi.f90
#   open_files.f90           read_arg.f90
#   read_conf.f90            send_det_mpi.f90         send_res_mpi.f90
#   term_loop.f90            prepare_it_oper.f90)
# target_link_libraries(bsr_breit4_mpi zconfls shells zcom)
endif()