include ${HOME}/DEFS_03

#  make  PROGRAM with name:
 
N  = test_int_pq
 
# Libraries to be searched:

L  =  $(DBS) $(ZCONFJJ) $(ZCOM) 

$(N): $(N).o 
	$(FC) -o $(N) $(link) $(N).o $(L) 
	rm *.o 

$(N).o: $(N).f90
	$(FC) $(comp) -fbacktrace -fcheck=bounds $(N).f90

