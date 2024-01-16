#Macro definition
FC = ifort -traceback
FFLAGS =  -O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
LDFLAGS =  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
           -Wl,--start-group \
           $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_core.a \
           $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
           -Wl,--end-group -lpthread -lm
OBJ = sub.o orient.o minpack.o geomtrans.o progdata.o diab.o hddata.o combinatorial.o \
			libutil.o libinternal.o localcoord.o io.o makesurf.o lm0.o
#end of Macro definition

surfgen.x: $(OBJ) surfgen0.o main.o
	$(FC) $(FFLAGS) $(OBJ) surfgen0.o main.o -o surfgen.x $(LDFLAGS) $(CUDAFLAGS)
clean:
	rm -f *.o *.mod a.out *.x *.a *.exe

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
