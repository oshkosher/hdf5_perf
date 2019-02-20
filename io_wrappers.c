#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdarg.h>
#include <pthread.h>
#include <mpi.h>

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
static int is_inside_mpi_file_open = 0;
static const char *mpi_file_open_filename = 0;

static int np = 0, rank = -1;

void io_wrappers_init();
void io_wrappers_report();
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
  
  is_init = 1;
}


void io_wrappers_report() {

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

  int result;

  /* This can be called before the shared library for MPI is available,
     and dlsym("MPI_Init") will fail, so avoid getting the pointers for all
     the functions in one initializer function. */

  
  if (needModeArg(flags)) {
    mode_t mode;
    va_list args;
    va_start(args, flags);
    mode = va_arg(args, mode_t);
    va_end(args);

    result = open_orig_fn(pathname, flags, mode);
    /*
    if (is_inside_mpi_file_open && result != -1) {
      fprintf(stderr, "[%d] open(\"%s\", 0%o, 0%o) -> fd=%d\n",
              rank, pathname, flags, mode, result);
    }
    */
  } else {
    result = open_orig_fn(pathname, flags);
    /*
    if (is_inside_mpi_file_open && result != -1) {
      fprintf(stderr, "[%d] open(\"%s\", 0%o) -> fd=%d\n",
              rank, pathname, flags, result);
    }
    */
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


int close(int fd) {
  static int (*close_orig_fn)(int fd) = 0;
  lookup((void**)&close_orig_fn, "close");

  if (isFileOpen(fd)) {
    fprintf(stderr, "[%d] close(%d)\n", rank, fd);
    setFileClosed(fd);
  }

  return close_orig_fn(fd);
}


ssize_t write(int fd, const void *buf, size_t count) {
  static ssize_t (*write_orig_fn)(int fd, const void *buf, size_t count) = 0;
  lookup((void**)&write_orig_fn, "write");

  ssize_t result = write_orig_fn(fd, buf, count);
  
  if (isFileOpen(fd)) {
    fprintf(stderr, "[%d] write(fd=%d, count=%lu)\n",
            rank, fd, (unsigned long)count);
    // printBacktrace();
  }

  return result;
}


ssize_t pwrite(int fd, const void *buf, size_t count, off_t offset) {
  static ssize_t (*pwrite_orig_fn)(int fd, const void *buf, size_t count, off_t offset) = 0;
  lookup((void**)&pwrite_orig_fn, "pwrite");
  double timer = MPI_Wtime();
  ssize_t result = pwrite_orig_fn(fd, buf, count, offset);
  timer = MPI_Wtime() - timer;
  
  if (isFileOpen(fd)) {
    fprintf(stderr, "[%d] pwrite(fd=%d, count=%lu, offset=%lu) %.6fs\n",
            rank, fd, (unsigned long)count, (unsigned long)offset, timer);
    /* printBacktrace(); */
  }

  return result;
}


ssize_t writev(int fd, const struct iovec *iov, int iovcnt) {
  static ssize_t (*writev_orig_fn)(int fd, const struct iovec *iov, int iovcnt) = 0;
  lookup((void**)&writev_orig_fn, "writev");

  ssize_t result = writev_orig_fn(fd, iov, iovcnt);
  
  if (isFileOpen(fd)) {
    fprintf(stderr, "[%d] writev(fd=%d, iovcnt=%d)\n", rank, fd, iovcnt);
    /* printBacktrace(); */
  }

  return result;
}


int MPI_Init(int *argc, char ***argv) {
  static int (*MPI_Init_orig_fn)(int *argc, char ***argv) = 0;
  int result;

  lookup((void**)&MPI_Init_orig_fn, "MPI_Init");

  result = MPI_Init_orig_fn(argc, argv);

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

  result = MPI_File_open_orig_fn(comm, filename, amode, info, fh);

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
  double timer = MPI_Wtime();

  result = fn(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);

  timer = MPI_Wtime() - timer;
  fprintf(stderr, "[%d] allgather %.6fs\n", rank, timer);

  return result;
}


int PMPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]) {
  static int (*fn)(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]) = 0;
  lookup((void**)&fn, "PMPI_Waitall");

  int result;
  double timer = MPI_Wtime();

  result = fn(count, array_of_requests, array_of_statuses);

  timer = MPI_Wtime() - timer;
  fprintf(stderr, "[%d] waitall %.6fs\n", rank, timer);

  return result;
}
                  

int PMPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request * request) {
  static int (*fn)(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request * request) = 0;
  lookup((void**)&fn, "PMPI_Irecv");

  int result;
  double timer = MPI_Wtime();

  result = fn(buf, count, datatype, source, tag, comm, request);

  timer = MPI_Wtime() - timer;
  fprintf(stderr, "[%d] irecv %.6fs\n", rank, timer);

  return result;
}


int PMPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
  static int (*fn)(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) = 0;
  lookup((void**)&fn, "PMPI_Send");
  int result;
  double timer = MPI_Wtime();

  result = fn(buf, count, datatype, dest, tag, comm);

  timer = MPI_Wtime() - timer;
  fprintf(stderr, "[%d] send %.6fs\n", rank, timer);

  return result;
}


int PMPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
  static int (*fn)(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) = 0;
  lookup((void**)&fn, "PMPI_Isend");
  int result;
  double timer = MPI_Wtime();

  result = fn(buf, count, datatype, dest, tag, comm, request);

  timer = MPI_Wtime() - timer;
  fprintf(stderr, "[%d] isend %.6fs\n", rank, timer);

  return result;
}

