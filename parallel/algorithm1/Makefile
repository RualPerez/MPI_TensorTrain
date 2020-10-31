# Sample Makefile
CC = mpicc
OPT = -g -O3
# WARN = -Wall
CFLAGS = $(OPT) # $(WARN) # the C compiler flags
OBJECTS = tt_contraction.o vandermonde.o init.o read_write.o read_write_int.o alloc.o operations.o

all: algorithm1 test_mse load_data

algorithm1: $(OBJECTS)
	$(CC) $(CFLAGS) -D ALG1 -o $@ $(OBJECTS) main.c

test_mse: $(OBJECTS)
	$(CC) $(CFLAGS) -D MSE -o $@ $(OBJECTS) main.c

load_data: $(OBJECTS)
	$(CC) $(CFLAGS) -D DATA -o $@ $(OBJECTS) main.c

clean:
	/bin/rm -f *.o core 

realclean: clean
	/bin/rm -f algorithm1 test_mse load_data

