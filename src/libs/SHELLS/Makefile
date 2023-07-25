include ${HOME}/DEFS_04

N  = libshells

S  =	mod_det_f3.f90           mod_det_f4.f90           mod_det_f5.f90           \
	detc_d2.f90              detc_d3.f90              detc_d4.f90              \
	detc_d5.f90              detc_d6.f90              detc_d7.f90              \
	detc_d8.f90              detc_f11.f90             detc_f12.f90             \
	detc_f3.f90              detc_f4.f90              detc_p2.f90              \
	detc_p3.f90              detc_p4.f90              detc_sh.f90              \
	detc_sh_LS.f90           det_d2.f90               det_d3.f90               \
	det_d4.f90               det_d5.f90               det_d6.f90               \
	det_d7.f90               det_d8.f90               det_f11.f90              \
	det_f12.f90              det_f2.f90                                        \
	det_p2.f90               det_p3.f90                                        \
	det_p4.f90               det_sh.f90               det_sh_LS.f90            \
	gen_f3.f90               gen_f4.f90               gen_f5.f90               \
	gen_f6.f90               gen_f7.f90               ML_id.f90                \
	zgen.f90

# det_f3.f90 det_f4.f90

C  =    det_exp.f

O = $(C:.f=.o) $(S:.f90=.o) 
 
$(N).a : $(O) 
	ar rv $(N).a $(O)
	mv $(N).a $(LIBDIR)
	rm *.o  
 
clean:
	rm -f *.o *~
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o .f
 
.f90.o:
	$(FC) $(compL) $<
 
.f.o:
	$(FC) $(compL) $<
 

