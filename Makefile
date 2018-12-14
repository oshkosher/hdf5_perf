default: all

EXECS=h5_collective
OPT=-g

# darshan-parser-nonzero /tmp/io.darshan | grep POSIX_ACCESS

# Blue Waters
ifeq "$(shell hostname | head -c 8)" "h2ologin"

LIB=
INC=
CC=cc $(OPT)
MPICC=cc $(OPT)
DARSHAN_RUNTIME=

else

HDF_HOME=/usr/local/hdf5-1.10.4
MPI_HOME=/usr/local/mpich-3.2.1
INC=-I$(HDF_HOME)/include -I$(MPI_HOME)/include
LIB=-L$(HDF_HOME)/lib -L$(MPI_HOME)/lib -lhdf5 -lmpi -Wl,-rpath -Wl,$(HDF_HOME)/lib
CC=gcc $(OPT) $(INC)
MPICC=mpicc $(OPT) $(INC)
# DARSHAN_LOGFILE=/tmp/io.darshan
# DARSHAN_RUNTIME=DXT_ENABLE_IO_TRACE=4 LD_PRELOAD=/usr/local/darshan3-runtime/lib/libdarshan.so DARSHAN_LOGFILE=$(DARSHAN_LOGFILE)

endif

all: $(EXECS)

h5_collective: h5_collective.c
	$(MPICC) $< $(LIB) -o $@

run: h5_collective
	mpirun -np 2 ./h5_collective /tmp/test_io 10 10
	rm -f /tmp/test_io

compare_write_sequence:
	grep 'X_POSIX[ 0]*write' meshio_large_file_details.txt | head -954 | cutfast -w -f 6,7 | ./sequence_extract.py busy

clean:
	rm -f $(EXECS) *.o *~
