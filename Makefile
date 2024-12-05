# librairies de SuiteSparse
L1 = SuiteSparse/UMFPACK/Lib/libumfpack.a
L2 = SuiteSparse/CHOLMOD/Lib/libcholmod.a 
L3 = SuiteSparse/AMD/Lib/libamd.a 
L4 = SuiteSparse/CAMD/Lib/libcamd.a  
L5 = SuiteSparse/COLAMD/Lib/libcolamd.a 
L6 = SuiteSparse/CCOLAMD/Lib/libccolamd.a 
L7 = SuiteSparse/metis-4.0/libmetis.a
L8 = SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a
LIB = $(L1) $(L2) $(L3) $(L4) $(L5) $(L6) $(L7) $(L8) -lm -lblas -llapack

COPT = -O3 -Wall -g

# list of agmg sources
listfiles = compute_density.o visualize.o prob.o time.o umfpack.o print_matrix.o compute_residue.o compute_charge.o utils.o

default: main

clean: 
	rm *.o 
	rm main

cleansrc:
	cd $(agmgdir); make clean

main: main.o $(listfiles)
	gcc -o main main.o $(listfiles) $(LIB)

main.o: main.c
	gcc $(COPT) -c main.c

umfpack.o: umfpack.c
	gcc $(COPT) -c $< -o  $@ -ISuiteSparse/UMFPACK/Include \
  -ISuiteSparse/SuiteSparse_config  -ISuiteSparse/AMD/Include

%.o: %.c
	gcc $(COPT) -c $< -o $@
