include ${HOME}/DEFS_04

 
N  = libzconfLS
 
S  =	MOD_param_LS.f90         MOD_symc.f90             MOD_symt.f90             \
	MOD_conf_LS.f90          MOD_orb_LS.f90           MOD_target.f90           \
	MOD_target_dc.f90        MOD_target_ion.f90       MOD_channel.f90          \
	MOD_channels.f90         MOD_channels_dc.f90      MOD_channels_ion.f90     \
	MOD_channel_dc.f90       MOD_overlaps_ion.f90     MOD_phys_orb.f90         \
	MOD_terms.f90            MOD_termexp.f90          MOD_spinorb.f90          \
	MOD_boef.f90             def_term.f90             idef_labelc.f90          \
	idef_nc.f90              ifind_channel.f90        incode_c.f90             \
	iort_conf.f90            iparity.f90              iterm_LS.f90             \
	label_c.f90              new_index.f90                                     \
	pri_conf.f90             r_closd.f90              r_confc.f90              \
	r_label.f90              r_term.f90                                        \
	test_c.f90               get_term_LS.f90                                   \
     incode_int.f90           idef_st.f90              r_expn.f90               \
     read_sol.f90             jopen.f90                SCT_tma_om.f90           \
     genterm.f90              sort_terms.f90           sort_states.f90          \
     read_cfgu.f90            prepare_it_oper.f90      prepare.f90              \
	idef_cme.f90             idet_fact.f90            add_it_oper.f90          \
	det_orbitals.f90         prepare_det_expn.f90


O = $(S:.f90=.o)
 
$(N).a : $(O)
	ar rv $(N).a $(O)
	mv $(N).a $(LIBDIR)
	rm *.o 
 
clean:
	rm -f *.o  *~ 
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(compL) $<
 
