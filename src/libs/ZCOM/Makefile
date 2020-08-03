include $(HOME)/DEFS_04

N  = libzcom
 
S  =	MOD_zconst.f90           MOD_coef.f90             MOD_det.f90   MOD_dets.f90    \
	MOD_id4.f90              MOD_int_lst.f90          MOD_rk.f90               \
	MOD_rk4.f90              MOD_tmat_ion.f90         MOD_tmat_list.f90        \
	MOD_tm_list.f90          MOD_tm_LS.f90            MOD_zoef.f90             \
	MOD_atoms_par.f90        MOD_atoms_LS.f90         MOD_orth.f90             \
	MOD_after.f90            MOD_new_dets.f90         MOD_new_defs.f90         \
	MOD_dip.f90              MOD_c_blocks.f90         MOD_orb_overlaps.f90     \
	angs_air.f90             ANG_dj_fact.f90          ANG_itra.f90             \
	ANG_kappa.f90            ANG_t_ls_jk.f90          ANG_zcklm.f90            \
	ANG_zclkl.f90            ANG_z_3j.f90             ANG_z_6j.f90             \
	ANG_z_9j.f90             check_file.f90           clean_a.f90              \
	conv_au.f90              det.f90                  full_mat.f90             \
	FUN_alegf.f90            FUN_bessel.f90           FUN_binom.f90            \
	FUN_cgamma.f90           FUN_dcoulfg.f90          FUN_dirac.f90            \
	FUN_gaussian.f90         FUN_hwf.f90              FUN_int_exp.f90          \
	FUN_legendre.f90                                  gauss_LA.f90             \
	gauss_points.f90         get_ij.f90               interv.f90               \
	inv.f90                  ipointer.f90                                      \
	LAPACK_calls.f90         num.f90                  read_arg.f90             \
	read_par.f90             recup_a.f90              rrtc.f90                 \
	SCT_coul_phase.f90       SCT_epahse.f90           SCT_fano_tma.f90         \
	SCT_iopen.f90            SCT_kma_om.f90           SCT_kma_tma.f90          \
	SCT_pri_mat.f90          SCT_sym_mat.f90          SCT_time_del.f90         \
	SCT_zkmat.f90            SCT_zrmat.f90            SCT_zsmat.f90            \
	sl_scale.f90             sorting_routines.f90     spectroscopic.f90        \
	splin3.f90               xlagr.f90                zis.f90                  \
	xAy.f90                  det_fact.f90             ndet_fact.f90      stop_mp.f90     \
	MOD_file.f90             FUN_zcoulfg90.f90        FUN_zcoulwf.f90    FUN_zcfg.f90    \
	chianti_gamma.f90        check_det.f90
	
 
O = $(S:.f90=.o)
 
$(N).a : $(O)
	ar rv $(N).a $(O)
	mv $(N).a  $(LIBDIR)
	rm *.o    
 
clean:
	rm -f *.o  *~
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o .f
 
.f90.o:
	$(FC) $(compL) $<
 
.f.o:
	$(FC) $(compL) $<
 
