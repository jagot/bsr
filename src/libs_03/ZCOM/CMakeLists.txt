set(ZCOM_SOURCES
  MOD_zconst.f90           MOD_coef_list.f90        MOD_det.f90   MOD_dets.f90
  MOD_id4.f90              MOD_int_lst.f90          MOD_rk.f90
  MOD_rk4.f90              MOD_tmat_ion.f90         MOD_tmat_list.f90
  MOD_tm_list.f90          MOD_tm_LS.f90            MOD_zoef.f90
  MOD_atoms_par.f90        MOD_atoms_LS.f90         MOD_orth.f90
  MOD_after.f90            MOD_new_dets.f90         MOD_new_defs.f90
  MOD_dip.f90              MOD_c_data.f90           MOD_obs.f90   MOD_file.f90
  angs_air.f90             ANG_dj_fact.f90          ANG_itra.f90
  ANG_kappa.f90            ANG_t_ls_jk.f90          ANG_zcklm.f90
  ANG_zclkl.f90            ANG_z_3j.f90             ANG_z_6j.f90
  ANG_z_9j.f90             check_file.f90           clean_a.f90
  conv_au.f90              det.f90                  full_mat.f90
  FUN_alegf.f90            FUN_bessel.f90           FUN_binom.f90
  FUN_cgamma.f90           FUN_dcoulfg.f90          FUN_dirac.f90
  FUN_gaussian.f90         FUN_hwf.f90              FUN_int_exp.f90
  FUN_legendre.f90         FUN_zcfg.f90             gauss_LA.f90
  gauss_points.f90         get_ij.f90               interv.f90
  inv.f90                  ipointer.f90
  LAPACK_calls.f90         num.f90                  read_arg.f90
  read_par.f90             recup_a.f90              rrtc.f90
  SCT_coul_phase.f90       SCT_epahse.f90           SCT_fano_tma.f90
  SCT_iopen.f90            SCT_kma_om.f90           SCT_kma_tma.f90
  SCT_pri_mat.f90          SCT_sym_mat.f90          SCT_time_del.f90
  SCT_zkmat.f90            SCT_zrmat.f90            SCT_zsmat.f90
  sl_scale.f90             sorting_routines.f90     spectroscopic.f90
  splin3.f90               xlagr.f90                zis.f90
  xAy.f90                  det_fact.f90             stop_mp.f90
  check_det_new.f90        det_fact_new.f90         index_tr.f90
  derivatives.f90)

# Gfortran no longer provides an implementation of iargc, so this
# simple routine is added when Gfortran is the compiler used.
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  list(APPEND ZCOM_SOURCES replace_iargc.f90)
endif()

add_library(zcom ${ZCOM_SOURCES})
