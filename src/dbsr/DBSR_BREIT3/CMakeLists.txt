add_executable(dbsr_breit3
  mod_dbsr_breit.f90   mod_nljm.f90         mod_term_exp.f90
  mod_coef.f90         mod_boef_list.f90
  dbsr_breit.f90       add_res.f90          conf_loop.f90
  det_exp.f90          det_me.f90           read_conf_jj.f90
  term_loop.f90        inf_dbsr_breit.f90   me_jj.f90)
target_link_libraries(dbsr_breit3 zconfjj zcom)

# # Conflict with utils?
# add_subdirectory(coef_jj)

add_subdirectory(det_jj)
