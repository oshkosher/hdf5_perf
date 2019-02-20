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
  if (fd >= open_files_len) return 0;
  return open_files[fd];
}


static void setFileClosed(int fd) {
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

void io_wrappers_init() {
  int fail = 0;

  if (is_init) return;
  /* fprintf(stderr, "io_wrappers_init tid = %x\n", (int) pthread_self()); */

  /*
  write_orig_fn = lookup("write", &fail);
  MPI_Init_orig_fn = lookup("MPI_Init", &fail);
  open_orig_fn = lookup("open", &fail);
  openat_orig_fn = lookup("openat", &fail);
  close_orig_fn = lookup("close", &fail);
  MPI_File_open_orig_fn = lookup("MPI_File_open", &fail);
  */
  if (fail) {
    fprintf(stderr, "ERROR: io_wrappers_init failed\n");
  } else {
    is_init = 1;
  }
}

void io_wrappers_report() {

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


static int getModeArg(int flags) {
  if (flags & O_CREAT) return 1;
#ifdef O_TMPFILE
  if (flags & O_TMPFILE) return 1;
#endif
  return 0;
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
  


int open(const char *pathname, int flags, ...) {
  int result;
  static int (*open_orig_fn)(const char *pathname, int flags, ...) = 0;

  /* This can be called before the shared library for MPI is available,
     and dlsym("MPI_Init") will fail, so avoid getting the pointers for all
     the functions in one initializer function. */

  lookup((void**)&open_orig_fn, "open");
  
  if (getModeArg(flags)) {
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

  if (getModeArg(flags)) {
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

  ssize_t result = pwrite_orig_fn(fd, buf, count, offset);
  
  if (isFileOpen(fd)) {
    fprintf(stderr, "[%d] pwrite(fd=%d, count=%lu, offset=%lu)\n",
            rank, fd, (unsigned long)count, (unsigned long)offset);
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



#if 0


int __wrap_puts(const char* str) {
  int result = __real_puts(str);
  printf("  call to puts with string of length %d\n", (int)strlen(str));
  //printBacktrace();
  fflush(stdout);
  return result;
}







ssize_t __wrap_write2(int fd, const void *buf, size_t count) {
  ssize_t result = __real_write(fd, buf, count);
  
  printf("  call to write, fd=%d, count=%d\n", fd, (int)count);

  return result;
}

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

