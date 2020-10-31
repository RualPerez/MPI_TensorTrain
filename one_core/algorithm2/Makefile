# Sample Makefile
CC = mpicc
OPT = -g -O3
WARN = -Wall
CFLAGS = $(OPT) $(WARN) # the C compiler flags
OBJECTS = alloc.o lapack_routines.o init.o kronecker_product.o operations.o read_write.o read_write_int.o reshape.o vandermonde.o

# lapack 
LIB_PATH = /zhome/0e/2/36189/teaching/02616/2019/libraries
LIBS = -L$(LIB_PATH)/lib -llapacke -llapack -lcblas -lblas -lgfortran -lm
CFLAGS += -I$(LIB_PATH)/include

algorithm2: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) main.c $(LIBS)

clean:
	/bin/rm -f *.o core 

realclean: clean
	/bin/rm -f algorithm2

