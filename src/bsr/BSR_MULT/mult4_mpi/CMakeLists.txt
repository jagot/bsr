add_executable(mult4_mpi
        mod_param.f90            mod_spinorb.f90          mod_termexp.f90
        mod_coef.f90             mod_zoef.f90             add_res.f90
        conf_loop.f90            def_ic.f90               det_mult.f90
        idet_fact.f90            mult.f90                 prepare.f90
        r_conf.f90               term_loop.f90
        stop_mpi.f90             read_arg.f90
        br_arg.f90               get_det_exp.f90          send_det_exp.f90
       conf_calc.f90            get_res_mpi.f90          send_res_mpi.f90)
target_link_libraries(mult4_mpi zconfls shells zcom MPI::MPI_Fortran)
