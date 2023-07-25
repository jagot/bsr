include ${HOME}/DEFS_03

# test RADIAL routines
test_rad: test_rad.f90
	$(FC) -o $@ $(link) $< $(RADIAL) $(ZCOM)
	mv $@ $(BINDIR)

