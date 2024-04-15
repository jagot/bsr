add_library(dbs
  DBS_nuclear.f90          DBS_grid.f90             DBS_gauss.f90
  DBS_orbitals_pq.f90      MOD_bk4.f90
  DBS_moments.f90          DBS_integrals.f90        DBS_core.f90   DBS_core_energy.f90
  DBS_dhl.f90              DBS_dhl_core.f90         DBS_debug.f90
  bdcwf.f90                bvalue2.f90              bav.f90
  bvalue_bm.f90            convert.f90              convol.f90
  def_Vnucl.f90            density.f90
  mrk_pppp.f90             mrk_pqpq.f90             mrk_qpqp.f90
  mrk_qqqq.f90             msk_ppqq.f90             msk_pqqp.f90
  pri_int.f90              quadr.f90                quadr_pq.f90
  read_pqbs.f90            rk_pq.f90                rk.f90
  sk.f90                   sk_pq.f90                sum_amb.f90
  tvm.f90                  update_hs.f90            zbsplvd.f90
  zint.f90                 zinty.f90                vinty.f90
  DBS_mult.f90             gen_hd.f90)
target_link_libraries(dbs zcom)
set_target_properties(dbs PROPERTIES Fortran_MODULE_DIRECTORY ${GLOBAL_Fortran_MODULE_DIRECTORY})
