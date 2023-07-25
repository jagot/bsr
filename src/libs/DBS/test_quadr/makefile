include ${HOME}/DEFS

#  make  PROGRAM with name:
 
N  = test_quadr
 
# Libraries to be searched:

L  =  $(DBS_new) $(ZCONFJJ) $(ZCOM) 

$(N): $(N).o 
	$(FC) -o $(N) $(link) $(N).o $(L) 
	rm *.o 

$(N).o: $(N).f90
	$(FC) $(comp) $(N).f90

