# Sample Makefile
CC = gcc
OPT = -g -O3
WARN = -Wall
CFLAGS = $(OPT) $(WARN) # the C compiler flags
OBJECTS = tt_contraction.o vandermonde.o init.o read_write.o read_write_int.o alloc.o

algorithm1: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) main.c

clean:
	/bin/rm -f *.o core 

realclean: clean
	/bin/rm -f algorithm1

