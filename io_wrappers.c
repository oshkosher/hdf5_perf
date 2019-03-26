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
#include <mpe_log.h>

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

#ifndef MPE_LOG_OK
#define MPE_Log_event(a,b,c)
#endif


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
static int is_inside_write = 0;
static const char *mpi_file_open_filename = 0;
static int is_inside_mpi_file_write = 0;
static int is_inside_mpi_file_read = 0;


static int np = 0, rank = -1;

/* put io counters (read/write/open/stat) first, then MPI calls, to quickly
   distinguish between acutal IO and communications overhead */
enum {
  IDX_WRITE,
  IDX_READ,

  IDX_ALLGATHER,
  IDX_ALLREDUCE,
  IDX_WAITALL,
  IDX_RECV,
  IDX_IRECV,
  IDX_SEND,
  IDX_ISEND,
  IDX_ALLTOALL,
  IDX_BARRIER,
  IDX_BCAST,

  IDX_COLL_WRITE,
  IDX_INDEP_WRITE,
  IDX_COLL_READ,
  IDX_INDEP_READ,

  IDX_END
};

#define IDX_FIRST_COMM IDX_ALLGATHER
#define IDX_FIRST_COLL_IO IDX_COLL_WRITE

const char *counter_names[] = {
  "write",
  "read",
  "allgather",
  "allreduce",
  "waitall",
  "recv",
  "irecv",
  "send",
  "isend",
  "alltoall",
  "barrier",
  "broadcast",
  "collective write",
  "independent write",
  "collective read",
  "independent read"
};

#define N_COUNTERS IDX_END

static int counters[N_COUNTERS];
static double timers[N_COUNTERS];
static uint64_t bytes[N_COUNTERS];

typedef struct {
  int counter_idx;
  int count;
  double timer;
} Counter;


static void io_wrappers_init(int do_text_log);
/* void io_wrappers_reset(); */ /* reset counters */
  /* collective routine; gather & print data */
/* void io_wrappers_report(FILE *outf); */
static void io_wrappers_end();  /* deallocate */
// static void printBacktrace();

// initialize logging and define states
static void init_mpe_logging();

/* log_ids stores MPE logging event ids. Entry N is the start event for
   state N/2, entry N+1 is the end event.
   These are set in init_mpe_logging() */
enum LogEventIds {
  STATE_COLL_WRITE = 0,
  STATE_WRITE = 2,
  STATE_WAITALL = 4,
  STATE_ALLTOALL = 6,
  STATE_ALLREDUCE = 8,
  STATE_ALLGATHER = 10,
  STATE_IRECV = 12,
  STATE_ISEND = 14,
  STATE_TYPE_COMMIT = 16,
  STATE_COLL_READ = 18,
  STATE_READ = 20,
  STATE_BARRIER = 22,
  STATE_BCAST = 24,
  STATE_INDEP_WRITE = 26,
  STATE_INDEP_READ = 28,
  
  LOG_EVENT_ID_END
};
int log_ids[LOG_EVENT_ID_END + 2];


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


void io_wrappers_init(int do_text_log) {
  if (is_init) return;
  
  assert(sizeof counter_names == sizeof(char*) * N_COUNTERS);
  assert(sizeof log_ids == sizeof(int) * (LOG_EVENT_ID_END + 2));

  t0 = 0;
  t0 = getTime();
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (do_text_log) {
    sprintf(log_filename, "io_wrappers.%d.out", rank);
    log = fopen(log_filename, "w");
    if (!log) {
      fprintf(stderr, "ERROR: failed to open %s for writing.\n", log_filename);
      log = stdout;
    } else {
      if (rank == 0) {
        printf("Writing io_wrappers.[0..%d].out\n", np-1);
      }
    }
  } else {
    log = NULL;
  }

#ifdef MPE_LOG_OK
  init_mpe_logging();
#endif

  is_init = 1;
}


/* reset counters */
void io_wrappers_reset() {
  memset(counters, 0, sizeof(int) * N_COUNTERS);
  memset(timers, 0, sizeof(double) * N_COUNTERS);
  /* memset(counters, 0, sizeof(Counter) * N_COUNTERS); */
  memset(bytes, 0, sizeof(uint64_t) * N_COUNTERS);
}


static int compareCounter(const void *av, const void *bv) {
  Counter *a = (Counter*) av;
  Counter *b = (Counter*) bv;

  if (a->count > 0) {
    if (b->count > 0) {
      /* both are nonzero, order by decreasing timer */
      if (a->timer < b->timer) {
        return 1;
      } else {
        return -1;
      }
    } else {
      /* only a is nonzero */
      return -1;
    }
  } else {
    if (b->count > 0) {
      /* only b is nonzero */
      return 1;
    } else {
      /* both are zero */
      return 0;
    }
  }
}


/* Fill 'ordered' with the values from counters[] and timers[], ordered
   by decreasing timer value. */
static void orderCounters(Counter ordered[N_COUNTERS],
                          const int counters[N_COUNTERS],
                          const double timers[N_COUNTERS]) {
  int i;
  /* copy counters and timers to one array of structures */
  for (i=0; i < N_COUNTERS; i++) {
    ordered[i].counter_idx = i;
    if (i < IDX_FIRST_COLL_IO) {
      ordered[i].count = counters[i];
      ordered[i].timer = timers[i];
    } else {
      /* don't count the time in collective io */
      ordered[i].count = 0;
      ordered[i].timer = 0;
    }
  }
  qsort(ordered, N_COUNTERS, sizeof(Counter), compareCounter);
}


/* collective routine; gather & print data */
void io_wrappers_report(FILE *outf) {
  int zero_counters[N_COUNTERS] = {0};
  double zero_timers[N_COUNTERS] = {0};

  int all_counters[N_COUNTERS];
  double all_timers[N_COUNTERS];
  uint64_t all_bytes[N_COUNTERS];

  /* print data for all processes */
  MPI_Reduce(counters, all_counters, N_COUNTERS, MPI_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(timers, all_timers, N_COUNTERS, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(bytes, all_bytes, N_COUNTERS, MPI_UINT64_T, MPI_SUM, 0,
             MPI_COMM_WORLD);

  if (rank == 0) {
    int idx = IDX_COLL_WRITE;
    int c = counters[idx];
    if (c) {
      fprintf(outf, "MPI_File_write %d calls, %.6fs, %ld bytes (raw write() %ld bytes)\n",
              c, timers[idx], (long)all_bytes[idx],
              (long)all_bytes[IDX_WRITE]);
    }

    idx = IDX_INDEP_WRITE;
    c = counters[idx];
    if (c) {
      fprintf(outf, "MPI_File_write %d calls, %.6fs, %ld bytes (raw write() %ld bytes)\n",
              c, timers[idx], (long)all_bytes[idx],
              (long)all_bytes[IDX_WRITE]);
    }

    idx = IDX_COLL_READ;
    c = counters[idx];
    if (c) {
      fprintf(outf, "MPI_File_read %d calls, %.6fs, %ld bytes (raw read() %ld bytes)\n",
              c, timers[idx], (long)all_bytes[idx],
              (long)all_bytes[IDX_READ]);
    }

    idx = IDX_INDEP_READ;
    c = counters[idx];
    if (c) {
      fprintf(outf, "MPI_File_read %d calls, %.6fs, %ld bytes (raw read() %ld bytes)\n",
              c, timers[idx], (long)all_bytes[idx],
              (long)all_bytes[IDX_READ]);
    }
  }

  int do_io, is_io;

  /* check whether this processes was actually doing I/O or just passing
     data to/from an IO processe */
  is_io = counters[IDX_WRITE] > 0 || counters[IDX_READ] > 0;

  // first print data from ranks that write, then ranks that don't
  for (int do_io=1; do_io >= 0; do_io--) {
    int do_me = is_io == do_io;
    int *counter_array = do_me ? counters : zero_counters;
    double *timer_array = do_me ? timers : zero_timers;
    int n_active;
    
    MPI_Reduce(&do_me, &n_active, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(counter_array, all_counters, N_COUNTERS, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(timer_array, all_timers, N_COUNTERS, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);

    double coll_io_time =
      all_timers[IDX_COLL_WRITE] + all_timers[IDX_COLL_READ] +
      all_timers[IDX_INDEP_WRITE] + all_timers[IDX_INDEP_READ];
    double raw_io_time = all_timers[IDX_WRITE] + all_timers[IDX_READ];

    if (rank == 0) {
      Counter ordered_counters[N_COUNTERS];
      orderCounters(ordered_counters, all_counters, all_timers);

      fprintf(outf, "%d %s processes, total time %.6fs", n_active,
              do_io ? "io" : "non-io", coll_io_time);

      if (do_io) {
        fprintf(outf, ", raw io %.6fs %.2f%%", 
                raw_io_time, 100.0 * raw_io_time / coll_io_time);
      }
      
      fprintf(outf, "\n");
      
      // double proc_sec = n_active * coll_io_time;

      Counter *c;
      double cumulative_pct = 0;
      for (c = ordered_counters; c != ordered_counters + N_COUNTERS; c++) {
        if (c->count) {
          double pct_time = 100.0 * c->timer / coll_io_time;
          cumulative_pct += pct_time;
          fprintf(outf, "  %s %.2f%% [%.2f%%] (%d call%s, %.6fs)\n",
                  counter_names[c->counter_idx], pct_time, cumulative_pct,
                  c->count, c->count==1 ? "" : "s", c->timer);
                  
        }
      }

    }
  }   
}


void io_wrappers_end() {
  if (!is_init) return;
  fclose(log);
}



#ifdef MPE_LOG_OK
// initialize logging and define states
static void init_mpe_logging() {
  int i;
  MPE_Init_log();

  for (i = 0; i < LOG_EVENT_ID_END; i += 2) {
    MPE_Log_get_state_eventIDs(log_ids + i, log_ids + i + 1);
  }

  if (rank == 0) {
    MPE_Describe_state(log_ids[STATE_COLL_WRITE], log_ids[STATE_COLL_WRITE+1],
                       "Collective write", "MidnightBlue");
    MPE_Describe_state(log_ids[STATE_INDEP_WRITE], log_ids[STATE_INDEP_WRITE+1],
                       "Independent write", "MidnightBlue");
    MPE_Describe_state(log_ids[STATE_COLL_READ], log_ids[STATE_COLL_READ+1],
                       "Collective read", "DarkGreen");
    MPE_Describe_state(log_ids[STATE_INDEP_READ], log_ids[STATE_INDEP_READ+1],
                       "Independent read", "DarkGreen");
    MPE_Describe_state(log_ids[STATE_WRITE], log_ids[STATE_WRITE+1],
                       "Write", "skyblue");
    MPE_Describe_state(log_ids[STATE_READ], log_ids[STATE_READ+1],
                       "Read", "green");
    MPE_Describe_state(log_ids[STATE_WAITALL], log_ids[STATE_WAITALL+1],
                       "Wait all", "white");
    MPE_Describe_state(log_ids[STATE_ALLTOALL], log_ids[STATE_ALLTOALL+1],
                       "All to all", "magenta");
    MPE_Describe_state(log_ids[STATE_ALLREDUCE], log_ids[STATE_ALLREDUCE+1],
                       "All reduce", "goldenrod");
    MPE_Describe_state(log_ids[STATE_ALLGATHER], log_ids[STATE_ALLGATHER+1],
                       "All gather", "magenta4");
    MPE_Describe_state(log_ids[STATE_BARRIER], log_ids[STATE_BARRIER+1],
                       "Barrier", "red");
    MPE_Describe_state(log_ids[STATE_BCAST], log_ids[STATE_BCAST+1],
                       "Broadcast", "wheat");
    MPE_Describe_state(log_ids[STATE_IRECV], log_ids[STATE_IRECV+1],
                       "Irecv", "gray25");
    MPE_Describe_state(log_ids[STATE_ISEND], log_ids[STATE_ISEND+1],
                       "Isend", "gray75");
    MPE_Describe_state(log_ids[STATE_TYPE_COMMIT],
                       log_ids[STATE_TYPE_COMMIT+1],
                       "Type commit", "ForestGreen");

  }
  
}

#endif



static int needModeArg(int flags) {
  if (flags & O_CREAT) return 1;
#ifdef O_TMPFILE
  if (flags & O_TMPFILE) return 1;
#endif
  return 0;
}

int open64(const char *pathname, int flags, ...) {
  static int (*fn)(const char *pathname, int flags, ...) = 0;
  lookup((void**)&fn, "open64");

  int fd;

  /* This can be called before the shared library for MPI is available,
     and dlsym("MPI_Init") will fail, so avoid getting the pointers for all
     the functions in one initializer function. */

  double start_time = getTime();
  mode_t mode = 0;
  if (needModeArg(flags)) {
    va_list args;
    va_start(args, flags);
    mode = va_arg(args, mode_t);
    va_end(args);

    fd = fn(pathname, flags, mode);
  } else {
    fd = fn(pathname, flags);
  }
  double end_time = getTime();

  if (log) fprintf(log, "open64(%s, %03o, %03o)->%d\n", pathname, flags, (int)mode, fd);
  
  // looks like this file was opened for MPI_File_open
  if (fd != -1
      // && !strcmp(pathname, mpi_file_open_filename)
      && is_inside_mpi_file_open) {
    double elapsed = end_time - start_time;
    if (log) fprintf(log, "%.6f open64 in MPI_File_open \"%s\" fd=%d %.6fs\n", getTime() - elapsed, pathname, fd, elapsed);
    setFileOpen(fd);
  }

  // fprintf(stderr, "open(%s)\n", pathname);

  return fd;
}


int open(const char *pathname, int flags, ...) {
  static int (*open_orig_fn)(const char *pathname, int flags, ...) = 0;
  lookup((void**)&open_orig_fn, "open");

  int fd;

  /* This can be called before the shared library for MPI is available,
     and dlsym("MPI_Init") will fail, so avoid getting the pointers for all
     the functions in one initializer function. */

  double start_time = getTime();
  mode_t mode = 0;
  if (needModeArg(flags)) {
    va_list args;
    va_start(args, flags);
    mode = va_arg(args, mode_t);
    va_end(args);

    fd = open_orig_fn(pathname, flags, mode);
  } else {
    fd = open_orig_fn(pathname, flags);
  }
  double end_time = getTime();

  if (log) fprintf(log, "open(%s, %03o, %03o)->%d\n", pathname, flags, (int)mode, fd);
  
  // looks like this file was opened for MPI_File_open
  if (fd != -1
      // && !strcmp(pathname, mpi_file_open_filename)
      && is_inside_mpi_file_open) {
    double elapsed = end_time - start_time;
    if (log) fprintf(log, "%.6f open in MPI_File_open \"%s\" fd=%d %.6fs\n", getTime() - elapsed, pathname, fd, elapsed);
    setFileOpen(fd);
  }

  // fprintf(stderr, "open(%s)\n", pathname);

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
    if (log) fprintf(log, "%.6f close(%d)\n", getTime(), fd);
    setFileClosed(fd);
  }

  return close_orig_fn(fd);
}


ssize_t write(int fd, const void *buf, size_t count) {
  static ssize_t (*fn)(int fd, const void *buf, size_t count) = 0;
  lookup((void**)&fn, "write");

  if (!isFileOpen(fd) || !is_inside_mpi_file_write) {

    return fn(fd, buf, count);

  } else {
  
    double start_time = getTime();
    MPE_Log_event(log_ids[STATE_WRITE], 0, NULL);
  
    ssize_t result = fn(fd, buf, count);

    MPE_Log_event(log_ids[STATE_WRITE+1], 0, NULL);
    double elapsed = getTime() - start_time;
    if (result > 0) bytes[IDX_WRITE] += result;

    int prev_is_inside = is_inside_write;
    is_inside_write = 1;
    /* fprintf(log, "%.6f write(fd=%d, count=%lu) %.6fs\n",
       start_time, fd, (unsigned long)count, elapsed); */
    counters[IDX_WRITE]++;
    timers[IDX_WRITE] += elapsed;
    is_inside_write = prev_is_inside;

    return result;

  }
}


ssize_t pwrite(int fd, const void *buf, size_t count, off_t offset) {
  static ssize_t (*fn)(int fd, const void *buf, size_t count, off_t offset) = 0;
  lookup((void**)&fn, "pwrite");
  
  if (!isFileOpen(fd) || !is_inside_mpi_file_write) {

    return fn(fd, buf, count, offset);

  } else {

    double start_time = getTime();
    MPE_Log_event(log_ids[STATE_WRITE], 0, NULL);

    ssize_t result = fn(fd, buf, count, offset);

    MPE_Log_event(log_ids[STATE_WRITE+1], 0, NULL);
    double elapsed = getTime() - start_time;
    if (result > 0) bytes[IDX_WRITE] += result;

    int prev_is_inside = is_inside_write;
    is_inside_write = 1;
    /* fprintf(log, "%.6f pwrite(fd=%d, count=%lu, offset=%lu) %.6fs\n",
            start_time, fd, (unsigned long)count, (unsigned long)offset,
            elapsed); */
    counters[IDX_WRITE]++;
    timers[IDX_WRITE] += elapsed;
    is_inside_write = prev_is_inside;
    /* printBacktrace(); */

    return result;

  }
}


ssize_t read(int fd, void *buf, size_t count) {
  static ssize_t (*fn)(int fd, void *buf, size_t count) = 0;
  lookup((void**)&fn, "read");
  
  if (!isFileOpen(fd) || !is_inside_mpi_file_read) {

    return fn(fd, buf, count);

  } else {
  
    double start_time = getTime();
    MPE_Log_event(log_ids[STATE_READ], 0, NULL);

    ssize_t result = fn(fd, buf, count);

    MPE_Log_event(log_ids[STATE_READ+1], 0, NULL);
    double elapsed = getTime() - start_time;
    if (result > 0) bytes[IDX_READ] += result;

    counters[IDX_READ]++;
    timers[IDX_READ] += elapsed;

    return result;

  }
}


ssize_t pread(int fd, void *buf, size_t count, off_t offset) {
  static ssize_t (*fn)(int fd, void *buf, size_t count, off_t offset) = 0;
  lookup((void**)&fn, "pread");
  
  if (!isFileOpen(fd) || !is_inside_mpi_file_read) {

    return fn(fd, buf, count, offset);

  } else {
  
    double start_time = getTime();
    MPE_Log_event(log_ids[STATE_READ], 0, NULL);

    ssize_t result = fn(fd, buf, count, offset);

    MPE_Log_event(log_ids[STATE_READ+1], 0, NULL);
    double elapsed = getTime() - start_time;
    if (result > 0) bytes[IDX_READ] += result;

    counters[IDX_READ]++;
    timers[IDX_READ] += elapsed;

    return result;

  }
}


int MPI_Init(int *argc, char ***argv) {
  static int (*fn)(int *argc, char ***argv) = 0;
  int result;

  lookup((void**)&fn, "MPI_Init");

  result = fn(argc, argv);
  /* fprintf(stderr, "wrapped MPI_Init\n"); */

  io_wrappers_init(0);
  
  return result;
}


int MPI_Finalize() {
  static int (*fn)(void) = 0;
  lookup((void**)&fn, "MPI_Finalize");

  if (open_files) {
    free(open_files);
    open_files = NULL;
    open_files_len = 0;
  }

  is_init = 0;

#ifdef MPE_LOG_OK
  const char *log_prefix = getenv("IO_WRAPPERS_LOG_PREFIX");
  if (!log_prefix)
    log_prefix = "io_wrappers";
  MPE_Finish_log(log_prefix);
  if (rank==0) printf("Wrote %s.clog2\n", log_prefix);
#endif
  
  int result = fn();

  return result;
}


int MPI_File_open(MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File *fh) {
  static int (*MPI_File_open_orig_fn)(MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File *fh) = 0;
  int result;

  lookup((void**)&MPI_File_open_orig_fn, "MPI_File_open");

  is_inside_mpi_file_open = 1;
  mpi_file_open_filename = filename;

  double start_time = getTime();
  if (log) fprintf(log, "%.6f MPI_File_open(%s) start\n", start_time, filename);

  result = MPI_File_open_orig_fn(comm, filename, amode, info, fh);

  double end_time = getTime();
  if (log) fprintf(log, "%.6f MPI_File_open(%s) end %.6f\n", end_time, filename, end_time - start_time);
          

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
  MPE_Log_event(log_ids[STATE_ALLGATHER], 0, NULL);

  result = fn(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);

  MPE_Log_event(log_ids[STATE_ALLGATHER+1], 0, NULL);
  double elapsed = getTime() - start_time;
  // fprintf(log, "%.6f allgather %.6fs\n", start_time, elapsed);
  counters[IDX_ALLGATHER]++;
  timers[IDX_ALLGATHER] += elapsed;

  return result;
}


int PMPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  static int (*fn)(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) = 0;
  lookup((void**)&fn, "PMPI_Allreduce");

  int result;
  double start_time = getTime();
  if (is_init)
    MPE_Log_event(log_ids[STATE_ALLREDUCE], 0, NULL);

  result = fn(sendbuf, recvbuf, count, datatype, op, comm);

  if (is_init)
    MPE_Log_event(log_ids[STATE_ALLREDUCE+1], 0, NULL);
  double elapsed = getTime() - start_time;
  // fprintf(stderr, "[%d] %.6f allreduce %.6fs\n", rank, start_time, elapsed);
  // fprintf(log, "%.6f allreduce %.6fs\n", start_time, elapsed);
  counters[IDX_ALLREDUCE]++;
  timers[IDX_ALLREDUCE] += elapsed;

  return result;
}



int PMPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]) {
  static int (*fn)(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]) = 0;
  lookup((void**)&fn, "PMPI_Waitall");

  int result;
  double start_time = getTime();
  MPE_Log_event(log_ids[STATE_WAITALL], 0, NULL);

  result = fn(count, array_of_requests, array_of_statuses);

  MPE_Log_event(log_ids[STATE_WAITALL+1], 0, NULL);
  double elapsed = getTime() - start_time;
  // fprintf(log, "%.6f waitall %.6fs\n", start_time, elapsed);
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
  // fprintf(log, "%.6f recv(count=%d, src=%d) %.6fs\n", start_time, count, source, elapsed);
  counters[IDX_RECV]++;
  timers[IDX_RECV] += elapsed;

  return result;
}
                  

int PMPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request * request) {
  static int (*fn)(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request * request) = 0;
  lookup((void**)&fn, "PMPI_Irecv");

  int result;
  double start_time = getTime();
  MPE_Log_event(log_ids[STATE_IRECV], 0, NULL);

  result = fn(buf, count, datatype, source, tag, comm, request);

  MPE_Log_event(log_ids[STATE_IRECV+1], 0, NULL);
  double elapsed = getTime() - start_time;
  // fprintf(log, "%.6f irecv(count=%d, src=%d, tag=%d) %.6fs\n", start_time, count, source, tag, elapsed);
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
  // fprintf(log, "%.6f send %.6fs\n", start_time, elapsed);
  counters[IDX_SEND]++;
  timers[IDX_SEND] += elapsed;

  return result;
}


int PMPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
  static int (*fn)(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) = 0;
  lookup((void**)&fn, "PMPI_Isend");
  int result;
  double start_time = getTime();
  MPE_Log_event(log_ids[STATE_ISEND], 0, NULL);

  result = fn(buf, count, datatype, dest, tag, comm, request);

  MPE_Log_event(log_ids[STATE_ISEND+1], 0, NULL);
  double elapsed = getTime() - start_time;
  // fprintf(log, "%.6f isend(count=%d, dest=%d, tag=%d) %.6fs\n", start_time, count, dest, tag, elapsed);
  counters[IDX_ISEND]++;
  timers[IDX_ISEND] += elapsed;

  return result;
}


int PMPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
  static int (*fn)(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) = 0;
  lookup((void**)&fn, "PMPI_Alltoall");
  int result;
  double start_time = getTime();
  MPE_Log_event(log_ids[STATE_ALLTOALL], 0, NULL);

  result = fn(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);

  MPE_Log_event(log_ids[STATE_ALLTOALL+1], 0, NULL);
  double elapsed = getTime() - start_time;
  // fprintf(log, "%.6f alltoall %.6fs\n", start_time, elapsed);
  counters[IDX_ALLTOALL]++;
  timers[IDX_ALLTOALL] += elapsed;

  return result;
}


int PMPI_Barrier(MPI_Comm comm) {
  static int (*fn)(MPI_Comm comm) = 0;
  lookup((void**)&fn, "PMPI_Barrier");

  int result;

  /* barrier seems to be called inside MPE_Finish_log, so logging inside
     that causes hangs */
  if (!is_init)
    return fn(comm);

  double start_time = getTime();
  MPE_Log_event(log_ids[STATE_BARRIER], 0, NULL);

  result = fn(comm);

  MPE_Log_event(log_ids[STATE_BARRIER+1], 0, NULL);
  double elapsed = getTime() - start_time;
  if (log) fprintf(log, "%.6f barrier %.6fs\n", start_time, elapsed);
  counters[IDX_BARRIER]++;
  timers[IDX_BARRIER] += elapsed;

  return result;
}


int PMPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
  static int (*fn)(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) = 0;
  lookup((void**)&fn, "PMPI_Bcast");

  int result;

  if (!is_init)
    return fn(buffer, count, datatype, root, comm);

  double start_time = getTime();
  MPE_Log_event(log_ids[STATE_BCAST], 0, NULL);

  result = fn(buffer, count, datatype, root, comm);

  MPE_Log_event(log_ids[STATE_BCAST+1], 0, NULL);
  double elapsed = getTime() - start_time;
  if (log) fprintf(log, "%.6f bcast %.6fs\n", start_time, elapsed);
  counters[IDX_BCAST]++;
  timers[IDX_BCAST] += elapsed;

  return result;
}



/* FYI, MPI_File is a pointer to a "struct ADIOI_FileD", defined in mpich-3.2.1/src/mpi/romio/adio/include/adio.h */
#ifdef __WINDOWS__
#define FDTYPE HANDLE
#else
#define FDTYPE int
#endif

#ifdef MPI_OFFSET_IS_INT
typedef int ADIO_Offset;
#else
typedef int64_t ADIO_Offset;
#endif


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

  /*
  if (rank == 0) {
    struct FileStruct *f = (struct FileStruct*) fh;
    fprintf(stderr, "MPI_File fh cookie=%x fd=%d d_mem=%u d_miniosz=%u blksize=%ld filename=%s\n", f->cookie, f->fd_sys, f->d_mem, f->d_miniosz, f->blksize, f->filename);
  }
  */

  int result;
  double start_time;
  uint64_t raw_io_bytes = bytes[IDX_WRITE];

  io_wrappers_reset();

  start_time = getTime();
  is_inside_mpi_file_write = 1;
  MPE_Log_event(log_ids[STATE_COLL_WRITE], 0, NULL);

  result = fn(fh, offset, buf, count, datatype, status);

  MPE_Log_event(log_ids[STATE_COLL_WRITE+1], 0, NULL);
  is_inside_mpi_file_write = 0;

  double elapsed = getTime() - start_time;
  counters[IDX_COLL_WRITE]++;
  timers[IDX_COLL_WRITE] += elapsed;
  raw_io_bytes = bytes[IDX_WRITE] - raw_io_bytes;

  /* compute the # of bytes written */
  int success_count;
  uint64_t mpi_io_bytes = 0;
  if (status != MPI_STATUS_IGNORE) {
    MPI_Get_count(status, datatype, &success_count);
  } else if (result == MPI_SUCCESS) {
    success_count = count;
  } else {
    success_count = 0;
  }
    
  if (success_count > 0) {
    MPI_Count dt_size;
    MPI_Type_size_x(datatype, &dt_size);
    mpi_io_bytes = success_count * dt_size;
    bytes[IDX_COLL_WRITE] += mpi_io_bytes;
  }

  if (log) fprintf(log, "%.6f file_write_at_all %ld bytes %.6fs write(%d calls, %ld bytes, %.6fs)\n", start_time, (long)mpi_io_bytes, elapsed, counters[IDX_WRITE], (long)raw_io_bytes, timers[IDX_WRITE]);

  if (rank == 0)
    printf("MPI_File_write_at_all report\n");
  io_wrappers_report(stdout);

  /*
  if (rank == 0) {
    printf("MPI_File_write_at_all total %.6fs=%.6f node*sec, io %.6fs (%.1f%%), "
           "comm %.6fs (%.1f%%)\n", elapsed, elapsed*np,
           io_time, elapsed == 0 ? 0 : 100.0 * io_time / (elapsed * np),
           comm_time, elapsed == 0 ? 0 : 100.0 * comm_time / (elapsed * np));
  }
  */

  return result;
}

  
int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void *buf, int count, MPI_Datatype datatype, MPI_Status *status) {
  static int (*fn)(MPI_File fh, MPI_Offset offset, const void *buf, int count, MPI_Datatype datatype, MPI_Status *status) = 0;
  lookup((void**)&fn, "MPI_File_write_at");

  /*
  if (rank == 0) {
    struct FileStruct *f = (struct FileStruct*) fh;
    fprintf(stderr, "MPI_File fh cookie=%x fd=%d d_mem=%u d_miniosz=%u blksize=%ld filename=%s\n", f->cookie, f->fd_sys, f->d_mem, f->d_miniosz, f->blksize, f->filename);
  }
  */

  int result;
  double start_time;
  uint64_t raw_io_bytes = bytes[IDX_WRITE];

  io_wrappers_reset();

  start_time = getTime();
  is_inside_mpi_file_write = 1;
  MPE_Log_event(log_ids[STATE_INDEP_WRITE], 0, NULL);

  result = fn(fh, offset, buf, count, datatype, status);

  MPE_Log_event(log_ids[STATE_INDEP_WRITE+1], 0, NULL);
  is_inside_mpi_file_write = 0;

  double elapsed = getTime() - start_time;
  counters[IDX_INDEP_WRITE]++;
  timers[IDX_INDEP_WRITE] += elapsed;
  raw_io_bytes = bytes[IDX_WRITE] - raw_io_bytes;

  /* compute the # of bytes written */
  int success_count;
  uint64_t mpi_io_bytes = 0;
  if (status != MPI_STATUS_IGNORE) {
    MPI_Get_count(status, datatype, &success_count);
  } else if (result == MPI_SUCCESS) {
    success_count = count;
  } else {
    success_count = 0;
  }
    
  if (success_count > 0) {
    MPI_Count dt_size;
    MPI_Type_size_x(datatype, &dt_size);
    mpi_io_bytes = success_count * dt_size;
    bytes[IDX_INDEP_WRITE] += mpi_io_bytes;
  }

  if (log) fprintf(log, "%.6f file_write_at %ld bytes %.6fs write(%d calls, %ld bytes, %.6fs)\n", start_time, (long)mpi_io_bytes, elapsed, counters[IDX_WRITE], (long)raw_io_bytes, timers[IDX_WRITE]);

  if (rank == 0)
    printf("MPI_File_write_at report\n");
  io_wrappers_report(stdout);

  /*
  if (rank == 0) {
    printf("MPI_File_write_at_all total %.6fs=%.6f node*sec, io %.6fs (%.1f%%), "
           "comm %.6fs (%.1f%%)\n", elapsed, elapsed*np,
           io_time, elapsed == 0 ? 0 : 100.0 * io_time / (elapsed * np),
           comm_time, elapsed == 0 ? 0 : 100.0 * comm_time / (elapsed * np));
  }
  */

  return result;
}


int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void * buf, int count,
                         MPI_Datatype datatype, MPI_Status *status) {
  static int (*fn)(MPI_File fh, MPI_Offset offset, void * buf, int count,
                   MPI_Datatype datatype, MPI_Status *status) = 0;
  lookup((void**)&fn, "MPI_File_read_at_all");

  int result;
  double start_time;
  uint64_t raw_io_bytes = bytes[IDX_READ];

  io_wrappers_reset();

  start_time = getTime();
  is_inside_mpi_file_read = 1;
  MPE_Log_event(log_ids[STATE_COLL_READ], 0, NULL);

  result = fn(fh, offset, buf, count, datatype, status);

  MPE_Log_event(log_ids[STATE_COLL_READ+1], 0, NULL);
  is_inside_mpi_file_read = 0;

  double elapsed = getTime() - start_time;
  counters[IDX_COLL_READ]++;
  timers[IDX_COLL_READ] += elapsed;
  raw_io_bytes = bytes[IDX_READ] - raw_io_bytes;

  /* compute the # of bytes read */
  int success_count;
  uint64_t mpi_io_bytes = 0;
  if (status != MPI_STATUS_IGNORE) {
    MPI_Get_count(status, datatype, &success_count);
  } else if (result == MPI_SUCCESS) {
    success_count = count;
  } else {
    success_count = 0;
  }
    
  if (success_count > 0) {
    MPI_Count dt_size;
    MPI_Type_size_x(datatype, &dt_size);
    mpi_io_bytes = success_count * dt_size;
    bytes[IDX_COLL_READ] += mpi_io_bytes;
  }

  if (log) fprintf(log, "%.6f file_read_at_all %ld bytes %.6fs read(%d calls, %ld bytes, %.6fs)\n", start_time, (long)mpi_io_bytes, elapsed, counters[IDX_READ], (long)raw_io_bytes, timers[IDX_READ]);

  if (rank == 0)
    printf("MPI_File_read_at_all report\n");
  io_wrappers_report(stdout);

  return result;
}  


int MPI_File_read_at(MPI_File fh, MPI_Offset offset, void * buf, int count,
                         MPI_Datatype datatype, MPI_Status *status) {
  static int (*fn)(MPI_File fh, MPI_Offset offset, void * buf, int count,
                   MPI_Datatype datatype, MPI_Status *status) = 0;
  lookup((void**)&fn, "MPI_File_read_at");

  int result;
  double start_time;
  uint64_t raw_io_bytes = bytes[IDX_READ];

  io_wrappers_reset();

  start_time = getTime();
  is_inside_mpi_file_read = 1;
  MPE_Log_event(log_ids[STATE_INDEP_READ], 0, NULL);

  result = fn(fh, offset, buf, count, datatype, status);

  MPE_Log_event(log_ids[STATE_INDEP_READ+1], 0, NULL);
  is_inside_mpi_file_read = 0;

  double elapsed = getTime() - start_time;
  counters[IDX_INDEP_READ]++;
  timers[IDX_INDEP_READ] += elapsed;
  raw_io_bytes = bytes[IDX_READ] - raw_io_bytes;

  /* compute the # of bytes read */
  int success_count;
  uint64_t mpi_io_bytes = 0;
  if (status != MPI_STATUS_IGNORE) {
    MPI_Get_count(status, datatype, &success_count);
  } else if (result == MPI_SUCCESS) {
    success_count = count;
  } else {
    success_count = 0;
  }
    
  if (success_count > 0) {
    MPI_Count dt_size;
    MPI_Type_size_x(datatype, &dt_size);
    mpi_io_bytes = success_count * dt_size;
    bytes[IDX_INDEP_READ] += mpi_io_bytes;
  }

  if (log) fprintf(log, "%.6f file_read_at %ld bytes %.6fs read(%d calls, %ld bytes, %.6fs)\n", start_time, (long)mpi_io_bytes, elapsed, counters[IDX_READ], (long)raw_io_bytes, timers[IDX_READ]);

  if (rank == 0)
    printf("MPI_File_read_at report\n");
  io_wrappers_report(stdout);

  return result;
}  


#if 0
int MPI_Type_commit(MPI_Datatype *datatype) {
  static int (*fn)(MPI_Datatype *datatype) = 0;
  lookup((void**)&fn, "MPI_Type_commit");
  
  MPE_Log_event(log_ids[STATE_TYPE_COMMIT], 0, NULL);

  int result = fn(datatype);

  MPE_Log_event(log_ids[STATE_TYPE_COMMIT+1], 0, NULL);

  return result;
}
#endif
