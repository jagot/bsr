include ${HOME}/DEFS_03

all:  print_bs  test_bs

# print B-splines in separate files for ploting
print_bs: print_bs.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

# test B-spline routines
test_bs: test_bs.f90
	$(FC) -o $@ $(link) $< $(BS) $(ZCOM) $(LAPACK)
	mv $@ $(BINDIR)

