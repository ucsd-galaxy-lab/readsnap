SHELL := /bin/bash

CC       =  mpicc
FC       =  $(CC)
HDF5INCL = -I/opt/hdf5/intel/mvapich2_ib/include -DH5_USE_16_API
HDF5LIB  = -L/opt/hdf5/intel/mvapich2_ib/lib -lhdf5 -lz
LIBS   = $(HDF5LIB)

CFLAGS = $(HDF5INCL)

EXEC = main
OBJS = main.c readsnap.c disk_finder.c calc_additional_parameters.c

.PHONY: help
help:
	cat Makefile

.PHONY: compile
compile: 
	cd ../ && make compile
	cp ../main readsnap

.PHONY: run
run: $(OBJS)
	source ./activate.sh && ibrun ./$(EXEC)

.PHONY: submit
submit:
	sbatch job_1-4
	sbatch job_2-4
	sbatch job_3-4
	sbatch job_4-4
	sbatch job_2-2
	sbatch job_4-2
	sbatch job_6-2
	sbatch job_1-8
	sbatch job_2-8


.PHONY: clean
clean:
	rm -f $(EXEC)
