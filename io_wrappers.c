#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdarg.h>
#include <inttypes.h>
#include <stdint.h>
#include <pthread.h>
#include <mpi.h>
#include "io_wrappers.h"

/*
  To wrap a function the uesr will call, wrap MPI_Xxx
  To wrap a function called inside MPI, wrap PMPI_Xxx
*/

#define __USE_GNU 1
#include <dlfcn.h>

#include <sys/uio.h>

/*
#define UNW_LOCAL_ONLY
#include <libunwind.h>
*/


/*
  -I$(UNWIND)/include
  -Wl,-wrap=write -Wl,-wrap=puts -Wl,-wrap=open
  -L$(UNWIND)/lib -lunwind
*/


// static int inside_wrapper = 0;

static int is_init = 0;
static FILE *log = 0;
static char log_filename[30];
static double t0 = 0;   // time 0, set in MPI_Init
static int is_inside_mpi_file_open = 0;
static const char *mpi_file_open_filename = 0;
static int is_inside_MPI_File_write_at_all = 0;


static int np = 0, rank = -1;

/* only count calls to write and pwrite inside MPI_File_write_at_all */
#if 0
static double timer_write = 0;
static int count_write = 0;

static double timer_Allgather = 0;
static int count_Allgather = 0;
static double timer_Waitall = 0;
static int count_Waitall = 0;
static double timer_Recv = 0;
static int count_Recv = 0;
static double timer_Irecv = 0;
static int count_Irecv = 0;
static double timer_Send = 0;
static int count_Send = 0;
static double timer_Isend = 0;
static int count_Isend = 0;
static double timer_Alltoall = 0;
static int count_Alltoall = 0;
#endif

/* put io counters (read/write/open/stat) first, then MPI calls */
enum {
  IDX_WRITE,

  IDX_ALLGATHER,
  IDX_WAITALL,
  IDX_RECV,
  IDX_IRECV,
  IDX_SEND,
  IDX_ISEND,
  IDX_ALLTOALL,

  IDX_END
};

#define IDX_FIRST_COMM IDX_ALLGATHER

const char *counter_names[] = {
  "write",
  "allgather",
  "waitall",
  "recv",
  "irecv",
  "send",
  "isend",
  "alltoall",
};

#define N_COUNTERS 8

static int counters[N_COUNTERS] = {0};
static double timers[N_COUNTERS] = {0};



void io_wrappers_init();
void io_wrappers_reset(); /* reset counters */
  /* collective routine; gather & print data */
void io_wrappers_report(FILE *outf, double *io_time, double *comm_time);
void io_wrappers_end();  /* deallocate */
// static void printBacktrace();

/* this will be an array indexed by file descriptor */
static int* open_files = NULL;
static int open_files_len = 0;

static void setFileOpen(int fd) {
  // check for unreasonable file descriptor
  assert(fd < 1000*1000*1000);

  // initialize
  if (!open_files) {
    open_files_len = 256;
    open_files = (int*) calloc(open_files_len, sizeof(int));
  }

  // grow
  if (fd >= open_files_len) {
    int *new_array, new_len;
    new_len = open_files_len;
    if (fd >= new_len)
      new_len = fd;
    new_array = (int*) calloc(new_len, sizeof(int));
    memcpy(new_array, open_files, sizeof(int) * open_files_len);
    open_files_len = new_len;
    free(open_files);
    open_files = new_array;
  }
  
  open_files[fd] = 1;
}


static int isFileOpen(int fd) {
  if (fd >= open_files_len || !open_files) return 0;
  return open_files[fd];
}


static void setFileClosed(int fd) {
  if (!open_files) return;
  if (fd < open_files_len)
    open_files[fd] = 0;
}


static void lookup(void **p, const char *name) {
  if (*p) return;
  *p = dlsym(RTLD_NEXT, name);
  if (!p) {
    fprintf(stderr, "ERROR: function \"%s\" not found by dlsym()\n", name);
  }
}


static double getTime() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (t.tv_sec + t.tv_usec * 0.000001) - t0;
}


#ifdef _UNWIND_H
static void printBacktrace() {
  unw_cursor_t cursor;
  unw_context_t context;

  // Initialize cursor to current frame for local unwinding.
  unw_getcontext(&context);
  unw_init_local(&cursor, &context);

  // Unwind frames one by one, going up the frame stack.
  while (unw_step(&cursor) > 0) {
    unw_word_t offset, pc;
    unw_get_reg(&cursor, UNW_REG_IP, &pc);
    if (pc == 0) {
      break;
    }
    printf("  0x%lx:", pc);

    char sym[256];
    if (unw_get_proc_name(&cursor, sym, sizeof(sym), &offset) == 0) {
      printf(" (%s+0x%lx)\n", sym, offset);
    } else {
      printf(" -- error: unable to obtain symbol name for this frame\n");
    }
  }
}
#endif


void io_wrappers_init() {
  if (is_init) return;
  
  assert(N_COUNTERS == IDX_END);
  assert(sizeof counter_names == sizeof(char*) * N_COUNTERS);

  is_init = 1;

  t0 = 0;
  t0 = getTime();
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  sprintf(log_filename, "io_wrappers.%d.out", rank);
  log = fopen(log_filename, "w");
  if (!log) {
    fprintf(stderr, "ERROR: failed to open %s for writing.\n", log_filename);
    log = stdout;
  } else {
    if (rank == 0) {
      fprintf(stderr, "writing io_wrappers.[0..%d].out\n", np-1);
    }
  }
}


/* reset counters */
void io_wrappers_reset() {
  /*
  timer_write = 0;
  count_write = 0;
  timer_Allgather = 0;
  count_Allgather = 0;
  timer_Waitall = 0;
  count_Waitall = 0;
  timer_Recv = 0;
  count_Recv = 0;
  timer_Irecv = 0;
  count_Irecv = 0;
  timer_Send = 0;
  count_Send = 0;
  timer_Isend = 0;
  count_Isend = 0;
  timer_Alltoall = 0;
  count_Alltoall = 0;
  */
  memset(counters, 0, sizeof(int) * N_COUNTERS);
  memset(timers, 0, sizeof(double) * N_COUNTERS);
}  


/* collective routine; gather & print data */
void io_wrappers_report(FILE *outf, double *io_time, double *comm_time) {
  int all_counters[N_COUNTERS];
  double all_timers[N_COUNTERS];

  MPI_Reduce(counters, all_counters, N_COUNTERS, MPI_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(timers, all_timers, N_COUNTERS, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);

  if (rank == 0) {
    int i;
    double total_time = 0;
    fprintf(outf, "cumulative counters\n");
    for (i=0; i < N_COUNTERS; i++) {
      if (all_counters[i]) {
        fprintf(outf, "  %s %d calls, %.6fs\n", counter_names[i],
                all_counters[i], all_timers[i]);
        total_time += all_timers[i];
      }
    }
    fprintf(outf, "comm overhead: %.6f\n", total_time);

    *io_time = *comm_time = 0;
  
    for (i=0; i < IDX_FIRST_COMM; i++)
      *io_time += all_timers[i];

    for (i=IDX_FIRST_COMM; i < N_COUNTERS; i++)
      *comm_time += all_timers[i];
  }
   
  /*
  if (count_write)
    fprintf(outf, "write() %d calls, %.6fs\n", count_write, timer_write);
  if (count_Allgather)
    fprintf(outf, "Allgather() %d calls, %.6fs\n", count_Allgather, timer_Allgather);
  if (count_Waitall)
    fprintf(outf, "Waitall() %d calls, %.6fs\n", count_Waitall, timer_Waitall);
  if (count_Recv)
    fprintf(outf, "Recv() %d calls, %.6fs\n", count_Recv, timer_Recv);
  if (count_Irecv)
    fprintf(outf, "Irecv() %d calls, %.6fs\n", count_Irecv, timer_Irecv);
  if (count_Send)
    fprintf(outf, "Send() %d calls, %.6fs\n", count_Send, timer_Send);
  if (count_Isend)
    fprintf(outf, "Isend() %d calls, %.6fs\n", count_Isend, timer_Isend);
  if (count_Alltoall)
    fprintf(outf, "Alltoall() %d calls, %.6fs\n", count_Alltoall, timer_Alltoall);
  
  fprintf(outf, "comm overhead: %.6f\n", timer_Allgather + timer_Waitall + timer_Recv + timer_Irecv + timer_Send + timer_Isend + timer_Alltoall);
  */
}


void io_wrappers_end() {
  if (!is_init) return;
  fclose(log);
}


static int needModeArg(int flags) {
  if (flags & O_CREAT) return 1;
#ifdef O_TMPFILE
  if (flags & O_TMPFILE) return 1;
#endif
  return 0;
}


int open(const char *pathname, int flags, ...) {
  static int (*open_orig_fn)(const char *pathname, int flags, ...) = 0;
  lookup((void**)&open_orig_fn, "open");

  int fd;

  /* This can be called before the shared library for MPI is available,
     and dlsym("MPI_Init") will fail, so avoid getting the pointers for all
     the functions in one initializer function. */

  double start_time = getTime();
  if (needModeArg(flags)) {
    mode_t mode;
    va_list args;
    va_start(args, flags);
    mode = va_arg(args, mode_t);
    va_end(args);

    fd = open_orig_fn(pathname, flags, mode);
  } else {
    fd = open_orig_fn(pathname, flags);
  }
  double end_time = getTime();
  
  // looks like this file was opened for MPI_File_open
  if (fd != -1
      && is_inside_mpi_file_open
      && !strcmp(pathname, mpi_file_open_filename)) {
    double elapsed = end_time - start_time;
    fprintf(log, "%.6f open in MPI_File_open \"%s\" fd=%d %.6fs\n",
            getTime() - elapsed, pathname, fd, elapsed);
    setFileOpen(fd);
  }

  return fd;
}


#if 0
int openat(int dirfd, const char *pathname, int flags, ...) {
  int result;
  static int (*openat_orig_fn)(int dirfd, const char *pathname, int flags, ...) = 0;
  lookup((void**)&openat_orig_fn, "openat");

  if (needModeArg(flags)) {
    mode_t mode;
    va_list args;
    va_start(args, flags);
    mode = va_arg(args, mode_t);
    va_end(args);

    result = openat_orig_fn(dirfd, pathname, flags, mode);
    if (is_inside_mpi_file_open && result != -1) {
      fprintf(stderr, "[%d] openat(\"%s\", 0%o, 0%o) -> fd=%d\n",
              rank, pathname, flags, mode, result);
    }
  } else {
    result = openat_orig_fn(dirfd, pathname, flags);
    if (is_inside_mpi_file_open && result != -1) {
      fprintf(stderr, "[%d] openat(\"%s\", 0%o) -> fd=%d\n",
              rank, pathname, flags, result);
    }
  }

  // looks like this file was opened for MPI_File_open
  if (result != -1
      && is_inside_mpi_file_open
      && !strcmp(pathname, mpi_file_open_filename)) {
    fprintf(stderr, "[%d] opened MPI_File \"%s\" fd=%d\n",
            rank, pathname, result);
    setFileOpen(result);
  }

  return result;
}
#endif


int close(int fd) {
  static int (*close_orig_fn)(int fd) = 0;
  lookup((void**)&close_orig_fn, "close");

  if (isFileOpen(fd)) {
    fprintf(log, "%.6f close(%d)\n", getTime(), fd);
    setFileClosed(fd);
  }

  return close_orig_fn(fd);
}


ssize_t write(int fd, const void *buf, size_t count) {
  static ssize_t (*write_orig_fn)(int fd, const void *buf, size_t count) = 0;
  lookup((void**)&write_orig_fn, "write");

  double start_time = getTime();
  ssize_t result = write_orig_fn(fd, buf, count);
  double elapsed = getTime() - start_time;
  
  if (isFileOpen(fd) && is_inside_MPI_File_write_at_all) {
    fprintf(log, "%.6f write(fd=%d, count=%lu) %.6fs\n",
            start_time, fd, (unsigned long)count, elapsed);
    counters[IDX_WRITE]++;
    timers[IDX_WRITE] += elapsed;
    /* printBacktrace(); */
  }

  return result;
}


ssize_t pwrite(int fd, const void *buf, size_t count, off_t offset) {
  static ssize_t (*pwrite_orig_fn)(int fd, const void *buf, size_t count, off_t offset) = 0;
  lookup((void**)&pwrite_orig_fn, "pwrite");
  double start_time = getTime();
  ssize_t result = pwrite_orig_fn(fd, buf, count, offset);
  double elapsed = getTime() - start_time;
  
  if (isFileOpen(fd) && is_inside_MPI_File_write_at_all) {
    fprintf(log, "%.6f pwrite(fd=%d, count=%lu, offset=%lu) %.6fs\n",
            start_time, fd, (unsigned long)count, (unsigned long)offset,
            elapsed);
    counters[IDX_WRITE]++;
    timers[IDX_WRITE] += elapsed;
    /* printBacktrace(); */
  }

  return result;
}


#if 0
ssize_t writev(int fd, const struct iovec *iov, int iovcnt) {
  static ssize_t (*writev_orig_fn)(int fd, const struct iovec *iov, int iovcnt) = 0;
  lookup((void**)&writev_orig_fn, "writev");

  ssize_t result = writev_orig_fn(fd, iov, iovcnt);
  
  if (isFileOpen(fd) && is_inside_MPI_File_write_at_all) {
    fprintf(stderr, "[%d] writev(fd=%d, iovcnt=%d)\n", rank, fd, iovcnt);
    /* printBacktrace(); */
  }

  return result;
}
#endif


int MPI_Init(int *argc, char ***argv) {
  static int (*MPI_Init_orig_fn)(int *argc, char ***argv) = 0;
  int result;

  lookup((void**)&MPI_Init_orig_fn, "MPI_Init");

  result = MPI_Init_orig_fn(argc, argv);

  io_wrappers_init();

  return result;
}


int MPI_Finalize() {
  static int (*MPI_Finalize_orig_fn)(void) = 0;
  lookup((void**)&MPI_Finalize_orig_fn, "MPI_Finalize");

  if (open_files) {
    free(open_files);
    open_files = NULL;
    open_files_len = 0;
  }

  return MPI_Finalize_orig_fn();
}


int MPI_File_open(MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File *fh) {
  static int (*MPI_File_open_orig_fn)(MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File *fh) = 0;
  int result;

  lookup((void**)&MPI_File_open_orig_fn, "MPI_File_open");

  is_inside_mpi_file_open = 1;
  mpi_file_open_filename = filename;

  double start_time = getTime();
  fprintf(log, "%.6f MPI_File_open(%s) start\n", start_time, filename);

  result = MPI_File_open_orig_fn(comm, filename, amode, info, fh);

  double end_time = getTime();
  fprintf(log, "%.6f MPI_File_open(%s) end %.6f\n", end_time, filename,
          end_time - start_time);

  is_inside_mpi_file_open = 0;
  mpi_file_open_filename = 0;

  return result;
}

/*
MPI_Allgather
MPI_Allreduce
MPI_Alltoall
MPI_Barrier
MPI_Bcast
MPI_Irecv
MPI_Isend
MPI_Recv
MPI_Send
MPI_Testall
MPI_Waitall
*/

int PMPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
  static int (*fn)(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) = 0;
  lookup((void**)&fn, "PMPI_Allgather");
  
  int result;
  double start_time = getTime();

  result = fn(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);

  double elapsed = getTime() - start_time;
  fprintf(log, "%.6f allgather %.6fs\n", start_time, elapsed);
  counters[IDX_ALLGATHER]++;
  timers[IDX_ALLGATHER] += elapsed;

  return result;
}


int PMPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]) {
  static int (*fn)(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]) = 0;
  lookup((void**)&fn, "PMPI_Waitall");

  int result;
  double start_time = getTime();

  result = fn(count, array_of_requests, array_of_statuses);

  double elapsed = getTime() - start_time;
  fprintf(log, "%.6f waitall %.6fs\n", start_time, elapsed);
  counters[IDX_WAITALL]++;
  timers[IDX_WAITALL] += elapsed;

  return result;
}


int PMPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, MPI_Status *status) {
  static int (*fn)(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) = 0;
  lookup((void**)&fn, "PMPI_Recv");

  int result;
  double start_time = getTime();

  result = fn(buf, count, datatype, source, tag, comm, status);

  double elapsed = getTime() - start_time;
  fprintf(log, "%.6f recv(count=%d, src=%d) %.6fs\n", start_time, count, source, elapsed);
  counters[IDX_RECV]++;
  timers[IDX_RECV] += elapsed;

  return result;
}
                  

int PMPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request * request) {
  static int (*fn)(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request * request) = 0;
  lookup((void**)&fn, "PMPI_Irecv");

  int result;
  double start_time = getTime();

  result = fn(buf, count, datatype, source, tag, comm, request);

  double elapsed = getTime() - start_time;
  fprintf(log, "%.6f irecv(count=%d, src=%d, tag=%d) %.6fs\n", start_time, count, source, tag, elapsed);
  counters[IDX_IRECV]++;
  timers[IDX_IRECV] += elapsed;

  return result;
}


int PMPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
  static int (*fn)(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) = 0;
  lookup((void**)&fn, "PMPI_Send");
  int result;
  double start_time = getTime();

  result = fn(buf, count, datatype, dest, tag, comm);

  double elapsed = getTime() - start_time;
  fprintf(log, "%.6f send %.6fs\n", start_time, elapsed);
  counters[IDX_SEND]++;
  timers[IDX_SEND] += elapsed;

  return result;
}


int PMPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
  static int (*fn)(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) = 0;
  lookup((void**)&fn, "PMPI_Isend");
  int result;
  double start_time = getTime();

  result = fn(buf, count, datatype, dest, tag, comm, request);

  double elapsed = getTime() - start_time;
  fprintf(log, "%.6f isend(count=%d, dest=%d, tag=%d) %.6fs\n", start_time, count, dest, tag, elapsed);
  counters[IDX_ISEND]++;
  timers[IDX_ISEND] += elapsed;

  return result;
}


int PMPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
  static int (*fn)(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) = 0;
  lookup((void**)&fn, "PMPI_Alltoall");
  int result;
  double start_time = getTime();

  result = fn(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);

  double elapsed = getTime() - start_time;
  fprintf(log, "%.6f alltoall %.6fs\n", start_time, elapsed);
  counters[IDX_ALLTOALL]++;
  timers[IDX_ALLTOALL] += elapsed;

  return result;
}


/* FYI, MPI_File is a pointer to a "struct ADIOI_FileD", defined in mpich-3.2.1/src/mpi/romio/adio/include/adio.h */
#ifdef __WINDOWS__
#define FDTYPE HANDLE
#else
#define FDTYPE int
#endif

typedef uint64_t ADIO_Offset;

struct FileStruct {
  int cookie;              /* for error checking */
  int fd_sys;              /* system file descriptor */
  int null_fd;          /* the null-device file descriptor: debug only (obviously)*/
  int fd_direct;           /* On XFS, this is used for direct I/O; 
                              fd_sys is used for buffered I/O */
  int direct_read;         /* flag; 1 means use direct read */
  int direct_write;        /* flag; 1 means use direct write  */
  /* direct I/O attributes */
  unsigned d_mem;          /* data buffer memory alignment */
  unsigned d_miniosz;      /* min xfer size, xfer size multiple,
                              and file seek offset alignment */
  long blksize;            /* some optimizations benefit from knowing
                              underlying block size */
  ADIO_Offset fp_ind;      /* individual file pointer in MPI-IO (in bytes)*/
  ADIO_Offset fp_sys_posn; /* current location of the system file-pointer
                              in bytes */
  void *fns;          /* struct of I/O functions to use */
  MPI_Comm comm;           /* communicator indicating who called open */
  int is_open;    /* deferred open: 0: not open yet 1: is open */
  int is_agg;              /* bool: if I am an aggregator */
  char *filename;          
  int file_system;         /* type of file system */
  int access_mode;         /* Access mode (sequential, append, etc.),
                              possibly modified to deal with
                              data sieving or deferred open*/
  int orig_access_mode;    /* Access mode provided by user: unmodified */
  ADIO_Offset disp;        /* reqd. for MPI-IO */
  MPI_Datatype etype;      /* reqd. for MPI-IO */
  MPI_Datatype filetype;   /* reqd. for MPI-IO */
  MPI_Count etype_size;          /* in bytes */
  void *hints;      /* structure containing fs-indep. info values */
  MPI_Info info;
};

  
int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf, int count, MPI_Datatype datatype, MPI_Status *status) {
  static int (*fn)(MPI_File fh, MPI_Offset offset, const void *buf, int count, MPI_Datatype datatype, MPI_Status *status) = 0;
  lookup((void**)&fn, "MPI_File_write_at_all");

  if (rank == 0) {
    struct FileStruct *f = (struct FileStruct*) fh;
    fprintf(stderr, "MPI_File fh cookie=%x fd=%d d_mem=%u d_miniosz=%u blksize=%ld filename=%s\n", f->cookie, f->fd_sys, f->d_mem, f->d_miniosz, f->blksize, f->filename);
  }

  int result;
  double start_time = getTime();

  io_wrappers_reset();
  is_inside_MPI_File_write_at_all = 1;

  result = fn(fh, offset, buf, count, datatype, status);

  is_inside_MPI_File_write_at_all = 0;

  double elapsed = getTime() - start_time;
  fprintf(log, "%.6f file_write_at_all %.6fs\n", start_time, elapsed);

  if (rank == 0)
    fprintf(stderr, "MPI_File_write_at_all report\n");
  double io_time, comm_time;
  io_wrappers_report(stdout, &io_time, &comm_time);

  if (rank == 0) {
    printf("MPI_File_write_at_all total %.6fs, io %.6fs (%.1f%%), "
           "comm %.6fs (%.1f%%)\n", elapsed,
           io_time, elapsed == 0 ? 0 : 100.0 * io_time / (elapsed * np),
           comm_time, elapsed == 0 ? 0 : 100.0 * comm_time / (elapsed * np));
  }

  return result;
}


