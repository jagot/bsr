add_library(radial
  mod_radial.f90           dyk.f90                  dzk.f90
  get_ecore.f90            grad.f90                 grad2.f90
  grad3.f90                hl.f90                   hlc.f90
  hlc_oo.f90               initr.f90                iwf.f90
  nk.f90                   p_deriv.f90              quadr.f90
  quads.f90                rk.f90                   rk_mass.f90
  rk_oo.f90                rlshft.f90               rw_wfn.f90
  tk.f90                   uk.f90                   vk.f90
  xk.f90                   zeta.f90                 zint.f90)
set_target_properties(radial PROPERTIES Fortran_MODULE_DIRECTORY ${GLOBAL_Fortran_MODULE_DIRECTORY})
