default: all

EXECS=h5_collective pc2orio2
OPT=-O3

# darshan-parser-nonzero /tmp/io.darshan | grep POSIX_ACCESS

# Blue Waters
ifeq "$(shell hostname | head -c 8)" "h2ologin"

HDF_HOME=/u/sciteam/karrels/Install/hdf5-1.10.4-cray
LIB=-L$(HDF_HOME)/lib -lhdf5 -Wl,-rpath -Wl,$(HDF_HOME)/lib
INC=-I$(HDF_HOME)/include
CC=cc $(INC) $(OPT)
MPICC=cc -DCRAY_MPI $(OPT) $(INC)
MPICXX=CC -DCRAY_MPI $(OPT) $(INC)
DARSHAN_RUNTIME=

else

HDF_HOME=/usr/local/hdf5
MPI_HOME=/usr/local/mpich
INC=-I$(HDF_HOME)/include
INC2=-I$(HDF_HOME)/include -I$(MPI_HOME)/include
LIB=-L$(HDF_HOME)/lib -lhdf5 -Wl,-rpath -Wl,$(HDF_HOME)/lib
LIB2=-L$(HDF_HOME)/lib -L$(MPI_HOME)/lib -lhdf5 -lmpi -Wl,-rpath -Wl,$(HDF_HOME)/lib
LIBDIRS=-L$(HDF_HOME)/lib -L$(MPI_HOME)/lib
CC=gcc $(OPT) $(INC)
MPICC=mpicc $(OPT) $(INC)
MPICXX=mpicxx $(OPT) $(INC)
# DARSHAN_LOGFILE=/tmp/io.darshan
# DARSHAN_RUNTIME=DXT_ENABLE_IO_TRACE=4 LD_PRELOAD=/usr/local/darshan3-runtime/lib/libdarshan.so DARSHAN_LOGFILE=$(DARSHAN_LOGFILE)

endif

all: $(EXECS)

h5_collective: h5_collective.c wrapper_fns.c
	$(MPICC) $< $(LIB) -o $@

UNWIND=/usr/local/libunwind-1.3
io_wrappers.o: io_wrappers.c
	gcc -c $< -I$(UNWIND)/include

pc2orio2: pc2orio2.cc
	$(MPICXX) $^ $(LIB) -lstdc++ -o $@

pc2orio2_wrapped: pc2orio2.cc io_wrappers.o
	g++ -static -g -pthread -I/usr/local/hdf5-1.10.4/include -I/usr/local/mpich-3.2.1/include pc2orio2.cc -L/usr/local/hdf5-1.10.4/lib -L/usr/local/mpich-3.2.1/lib -lhdf5 -lmpi -Wl,-wrap=write -Wl,-wrap=pwrite -Wl,-wrap=writev -Wl,-wrap=puts -Wl,-wrap=open -Wl,-wrap=close -L/usr/local/libunwind-1.3/lib io_wrappers.o -lunwind -lz -ldl -lstdc++ -lrt -o pc2orio2


# TAU_OPT=-tau_options=-optTauSelectFile=pc2orio2.inst
# TAU_OPT=-tau_options=-optCompInst
# TAU_OPT=-tau_options=–optTrackIO
pc2orio2_tau: pc2orio2.cc
	rm -f pc2orio2.o
	tau_cxx.sh $(TAU_OPT) -c $< $(INC)
	tau_cxx.sh pc2orio2.o -o $@ $(LIBDIRS) -lhdf5 -lmpi

# g++ -o $@ pc2orio2.o $(LIBDIRS) -lhdf5  `tau_cc.sh -tau:showlibs`

# /usr/bin/g++ -L/usr/local/hdf5-1.10.4/lib -L/usr/local/mpich-3.2.1/lib -L/usr/local/tau-2.28/x86_64/lib -L/usr/local/tau-2.28/x86_64/binutils-2.23.2/lib -L/usr/local/tau-2.28/x86_64/binutils-2.23.2/lib64 -L/usr/local/tau-2.28/x86_64/libunwind-1.3-rc1-gcc/lib pc2orio2.o -lhdf5 -lTauMpi-mpi-pdt-trace -ltau-mpi-pdt-trace -lmpi -lbfd -liberty -lz -ldl -lrt -lunwind -lm -lstdc++ -lgcc_s -g -Wl,--export-dynamic -o pc2orio2_tau

# tau_cxx.sh $< -o $@ $(INC) $(LIBDIRS) -lhdf5 -lmpi

# -optLinking="-lhdf5 -lmpi" 

# /usr/bin/g++ -L/usr/local/hdf5-1.10.4/lib -L/usr/local/mpich-3.2.1/lib -L/usr/local/tau-2.28/x86_64/lib -L/usr/local/tau-2.28/x86_64/binutils-2.23.2/lib -L/usr/local/tau-2.28/x86_64/binutils-2.23.2/lib64 -L/usr/local/tau-2.28/x86_64/libunwind-1.3-rc1-gcc/lib pc2orio2.o -lhdf5 -lTauMpi-mpi-pdt-trace -ltau-mpi-pdt-trace -lmpi -lbfd -liberty -lz -ldl -lrt -lunwind -lm -lstdc++ -lgcc_s -g -Wl,--export-dynamic -o pc2orio2_tau

run: pc2orio2
	mpirun -np 2 ./pc2orio2 /tmp/test_io 100,100,100 2,1,1 3
	rm -f /tmp/test_io*

tau: pc2orio2_tau
	rm -f *.edf *.trc
	mpirun -np 2 ./pc2orio2_tau /tmp/test_io 100,100,100 2,1,1 3
	tau_treemerge.pl
	tau_convert -dump tau.trc tau.edf | head -100

clean:
	rm -f $(EXECS) *.o *~ *.trc *.edf *.slog2
