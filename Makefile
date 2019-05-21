SHELL := /bin/bash

CC       =  mpicc
FC       =  $(CC)
HDF5INCL = -I/opt/hdf5/intel/mvapich2_ib/include -DH5_USE_16_API
HDF5LIB  = -L/opt/hdf5/intel/mvapich2_ib/lib -lhdf5 -lz
LIBS   = $(HDF5LIB)

CFLAGS = $(HDF5INCL)

EXEC = readsnap
OBJS = readsnap.c

.PHONY: help
help:
	cat Makefile

.PHONY: compile
compile: $(OBJS)
	source ./activate.sh && $(FC) $(CFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

.PHONY: run
run: $(OBJS)
	source ./activate.sh && ibrun ./$(EXEC)

.PHONY: clean
clean:
	rm -f $(EXEC)
