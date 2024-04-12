add_executable(jj2ls
  mod_jj2ls.f90            lq_jj_ls.f90             check_jj_shells.f90
  check_ovl.f90            cm2j.f90                 get_LS_conf.f90
  get_LS_terms.f90         get_tr_coef.f90          g_factor.f90
  jj2ls.f90                make_couplings.f90       transform_jfile.f90
  recup_3j.f90             recup_jj_shells.f90      sub_nlLS.f90)
target_link_libraries(jj2ls zconfjj zconfls zcom)
