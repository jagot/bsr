include ${HOME}/DEFS

#  make  PROGRAM with name:
 
N  = test_grid_DBS_new
 
# Libraries to be searched:

L  =  $(DBS_new) $(ZCONFJJ) $(ZCOM) 

$(N): $(N).o 
	$(FC) -o $(N) $(link) $(N).o $(L) 
	mv $(N) $(BINDIR)
	rm *.o 

$(N).o: $(N).f90
	$(FC) $(comp) $(N).f90

