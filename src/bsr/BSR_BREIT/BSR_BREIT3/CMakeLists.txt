add_executable(bsr_breit3
  mod_param.f90            mod_termexp.f90          mod_spinorb.f90
  mod_boef.f90             mod_coef.f90             mod_zoef.f90
  add_it_oper.f90          add_res.f90              bsr_breit.f90
  conf_loop.f90            det_btreit.f90           idef_cme.f90
  idet_fact.f90            open_br.f90              prepare.f90
  read_arg.f90             read_conf.f90            term_loop.f90
  zno_0ee.f90              zno_1ee.f90              zno_2ee.f90
  zno_breit.f90            iort_conf.f90)
target_link_libraries(bsr_breit3 zconfls shells zcom)

add_subdirectory(coef_tab3)
