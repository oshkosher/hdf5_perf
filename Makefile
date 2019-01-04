default: all

EXECS=h5_collective pc2orio2
OPT=-O3

# darshan-parser-nonzero /tmp/io.darshan | grep POSIX_ACCESS

# Blue Waters
ifeq "$(shell hostname | head -c 8)" "h2ologin"

LIB=
INC=
CC=cc $(OPT)
MPICC=cc $(OPT)
MPICXX=CC $(OPT)
DARSHAN_RUNTIME=

else

HDF_HOME=/usr/local/hdf5-1.10.4
MPI_HOME=/usr/local/mpich-3.2.1
INC=-I$(HDF_HOME)/include -I$(MPI_HOME)/include
LIB=-L$(HDF_HOME)/lib -L$(MPI_HOME)/lib -lhdf5 -lmpi -Wl,-rpath -Wl,$(HDF_HOME)/lib
CC=gcc $(OPT) $(INC)
MPICC=mpicc $(OPT) $(INC)
MPICXX=mpicxx $(OPT) $(INC)
# DARSHAN_LOGFILE=/tmp/io.darshan
# DARSHAN_RUNTIME=DXT_ENABLE_IO_TRACE=4 LD_PRELOAD=/usr/local/darshan3-runtime/lib/libdarshan.so DARSHAN_LOGFILE=$(DARSHAN_LOGFILE)

endif

all: $(EXECS)

h5_collective: h5_collective.c wrapper_fns.c
	$(MPICC) $< $(LIB) -o $@

pc2orio2: pc2orio2.cc
	$(MPICC) $< $(LIB) -lstdc++ -o $@

run: pc2orio2
	mpirun -np 2 ./pc2orio2 /tmp/test_io 100,100,100 2,1,1 3
	rm -f /tmp/test_io*

clean:
	rm -f $(EXECS) *.o *~
