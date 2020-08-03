add_library(bs
  spline_param.f90         spline_atomic.f90        spline_grid.f90
  spline_galerkin.f90      spline_orbitals.f90      spline_densities.f90
  spline_integrals.f90     spline_moments.f90       spline_hl.f90
  spline_rk_integrals.f90  spline_slater.f90        azl.f90
  bav.f90                  bcore.f90                bhl.f90
  bhwf.f90                 bs_cwf.f90               bvalue2.f90
  bvmv.f90                 bxv.f90                  bzk.f90
  convol.f90               coulom.f90               define_grid.f90
  define_spline.f90        density.f90
  dinty.f90                dvx.f90                  facdyk.f90
  facdzk.f90               facsb.f90                facsbl.f90
  gauss.f90                grad.f90                 hlc.f90
  hlm.f90                  int_de.f90               int_v.f90
  minty.f90                mk.f90                   mkc.f90
  mkgrid2.f90              mky.f90                  mk_moments.f90
  mmk_cell.f90             mmk_diff.f90             mnk_cell.f90
  mnk_diff.f90             moments.f90              mqk_cell.f90
  mrk_cell.f90             mrk_diff.f90             mrm.f90
  mtk_cell.f90             mtk_diff.f90             mvc.f90
  mvcv.f90                 mvk_cell.f90             mvk_diff.f90
  mwk_cell.f90             nk.f90                   nkc.f90
  nky.f90                  nk_moments.f90           qk.f90
  qkc.f90                  qky.f90                  qk_moments.f90
  quadr.f90                rk.f90                   rkc.f90
  rky.f90                  rk_moments.f90           r_bwfn.f90
  splin3.f90               sum_amb.f90              tk.f90
  tkc.f90                  tky.f90                  tk_moments.f90
  vbsplvd.f90              vinty.f90                vk.f90
  vkc.f90                  vky.f90                  vk_moments.f90
  wk.f90                   wky.f90                  wk_moments.f90
  ykf.f90                  yval.f90                 zeta_y.f90
  zinty.f90                memory_splines.f90)
