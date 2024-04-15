add_executable(coef_tab3 mod_param.f90   coef_tab3.f90   pri_coef.f90  rr_conf.f90  IDET_SIMP1.f90)
target_link_libraries(coef_tab3 zconfls zcom)
