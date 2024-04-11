# extract seperate states from "bound.nnn" file:
add_executable(bound_bsw bound_bsw.f90)
target_link_libraries(bound_bsw bs zconfls zcom ${BLAS_LAPACK_LIBRARIES})

# list of states in the "bound.nnn" files:
add_executable(bound_tab bound_tab.f90)
target_link_libraries(bound_tab zcom)

# hydrogen-like bound and pseudocontinuum states for given B-soline basis:
add_executable(bs_hwf bs_hwf.f90)
target_link_libraries(bs_hwf bs zcom ${BLAS_LAPACK_LIBRARIES})

# calculation of dipole matrix elements between R-matrx states saved in rsol.nnn files
add_executable(bsr_dd bsr_dd.f90)
target_link_libraries(bsr_dd zconfls zcom)

# merging  1.bsw + 2.bsw + 3.bsw + ... --> res.bsw
add_executable(bsw123 bsw123.f90)
target_link_libraries(bsw123 bs zcom ${BLAS_LAPACK_LIBRARIES})

# connverts the bsw-files with orbitals in B-spline repersentation 
# into tab-files suitable for grafic display
add_executable(bsw_tab bsw_tab.f90)
target_link_libraries(bsw_tab bs zcom ${BLAS_LAPACK_LIBRARIES})

# connvert the B-spline results in the MCHF format
add_executable(bsw_w bsw_w.f90)
target_link_libraries(bsw_w bs zcom ${BLAS_LAPACK_LIBRARIES})

# target, rsol.nnn  --> rsol.all
add_executable(collect_rsol collect_rsol.f90)
target_link_libraries(collect_rsol zcom)

# # target, threshoulds -> target.exp
# add_executable(convert_target convert_target.f90)
# target_link_libraries(convert_target zconfls zcom)

# # H.DAT  -->  f-values
# add_executable(f_values f_values.f90)
# target_link_libraries(f_values zconfls zcom)

# different photoionization paramters
add_executable(photo_collect photo_collect.f90)
target_link_libraries(photo_collect zcom)

# # different photoionization paramters
# add_executable(photo_tab_new photo_tab_new.f90)
# target_link_libraries(photo_tab_new zconfls zcom)

# prints B-splines in separate files:
add_executable(print_bs print_bs.f90)
target_link_libraries(print_bs bs zcom ${BLAS_LAPACK_LIBRARIES})

# # H.DAT --> dipole matrix elements
# add_executable(s_values s_values.f90)
# target_link_libraries(s_values zconfls zcom)

# target states information
add_executable(target_states target_states.f90)
target_link_libraries(target_states zconfls zcom)

# Generate B-spline representation for orbitals 
# given in MCHF package format, w-files 
add_executable(w_bsw w_bsw.f90)
target_link_libraries(w_bsw bs zcom ${BLAS_LAPACK_LIBRARIES})

# Make H.DAT from h.nnn file
add_executable(sum_hh_jK sum_hh.f90)
target_link_libraries(sum_hh_jK zcom)

# Make H.DAT from h.nnn file for jj-coupling case
add_executable(sum_hh_jj sum_hh_jj.f90)
target_link_libraries(sum_hh_jj zcom)

# Make top-up executable case
add_executable(sec_top_TM sec_top_TM.f90)
target_link_libraries(sec_top_TM zconfls zcom)

# Make top-up executable case no "TM"
add_executable(sec_top sec_top.f90)
target_link_libraries(sec_top zconfls zcom)

# Convert tma to tmb
add_executable(tma_tmb tma_tmb.f90)
target_link_libraries(tma_tmb zconfls zcom)

# Convert tmb to omb
add_executable(tmb_omb tmb_omb.f90)
target_link_libraries(tmb_omb zconfls zcom)

# Another top-up executable 
add_executable(sec_omn sec_omn.f90)
target_link_libraries(sec_omn zconfls zcom)

# Make .tma and .om from .kma
add_executable(kma_om kma_om.f90)
target_link_libraries(kma_om zconfls zcom)

# Make .omb from .oma
add_executable(oma_omb oma_omb.f90)
target_link_libraries(oma_omb zconfls zcom)

# Make f-value file
add_executable(f_values_jj f_values_jj.f90)
target_link_libraries(f_values_jj zconfjj zcom)

#  f-value calculations between states in bound.nnn files
add_executable(zf_cc_bsr zf_cc_bsr.f90)
target_link_libraries(zf_cc_bsr zconfls zcom)

#  f-value calculations between set of states in the c- and bsw-files
add_executable(zf_bb_bsr zf_bb_bsr.f90)
target_link_libraries(zf_bb_bsr zconfls zcom)

add_executable(target_kpert target_kpert.f90)
target_link_libraries(target_kpert zconfls zcom)

add_executable(add_farm add_farm.f90)
target_link_libraries(add_farm zconfls zcom)

add_executable(add_stgf add_stgf.f90)
target_link_libraries(add_stgf zconfls zcom)

add_executable(add_stgfb add_stgfb.f90)
target_link_libraries(add_stgfb zconfls zcom)

add_executable(sec_om sec_om.f90)
target_link_libraries(sec_om zconfls zcom)

add_executable(ak_tab_jj ak_tab_jj.f90)
target_link_libraries(ak_tab_jj zconfjj zcom)

add_executable(cm_j cm_j.f90)
target_link_libraries(cm_j zconfjj zcom)

add_executable(genjconf genjconf.f90)
target_link_libraries(genjconf zconfjj zcom)

add_executable(genjterm genjterm.f90)
target_link_libraries(genjterm zconfjj zcom)

add_executable(jcfile jcfile.f90)
target_link_libraries(jcfile zconfjj zcom)

add_executable(jfile jfile.f90)
target_link_libraries(jfile zconfjj zcom)

add_executable(jj_list jj_list.f90)
target_link_libraries(jj_list zconfjj zcom)

add_executable(mix_c mix_c.f90)
target_link_libraries(mix_c zconfjj zcom)

add_executable(order_jj order_jj.f90)
target_link_libraries(order_jj zconfjj zcom)
