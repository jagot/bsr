message(WARNING "bsr_breit4 does not build correctly; requires version 4 libraries to build")
# add_executable(bsr_breit4
#   MOD_bsr_breit.f90                                              
#   BSR_BREIT.f90            breit_me.f90             add_res.f90              
#   conf_loop.f90            open_files.f90           read_arg.f90             
#   read_conf.f90            term_loop.f90            prepare_it_oper.f90)
# target_link_libraries(bsr_breit4 zconfls4 shells4 zcom4)

add_subdirectory(coef_tab4)
