#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/uio.h>

#define UNW_LOCAL_ONLY
#include <libunwind.h>

/*
  -I$(UNWIND)/include
  -Wl,-wrap=write -Wl,-wrap=puts -Wl,-wrap=open
  -L$(UNWIND)/lib -lunwind
*/


static int inside_wrapper = 0;

ssize_t __real_write(int fd, const void *buf, size_t count);
ssize_t __real_pwrite(int fd, const void *buf, size_t count, off_t offset);
ssize_t __real_writev(int fd, const struct iovec *iov, int iovcnt);
int __real_puts(const char* str);
int __real_open(const char *pathname, int flags, mode_t mode);
int __real_close(int fd);

static void printBacktrace();

int __wrap_open(const char *pathname, int flags, mode_t mode) {
  int prev_inside_wrapper = inside_wrapper;
  inside_wrapper = 1;
  int result = __real_open(pathname, flags, mode);
  if (!prev_inside_wrapper) {
    printf("  open(\"%s\", 0%o, 0%o) -> fd=%d\n", pathname, flags, mode, result);
    /* printBacktrace(); */
  }
  inside_wrapper = prev_inside_wrapper;
  return result;
}


int __wrap_close(int fd) {
  int prev_inside_wrapper = inside_wrapper;
  inside_wrapper = 1;
  if (!prev_inside_wrapper) {
    printf("  close(%d)\n", fd);
  }
  inside_wrapper = prev_inside_wrapper;
  return __real_close(fd);
}



int __wrap_puts(const char* str) {
  int result = __real_puts(str);
  printf("  call to puts with string of length %d\n", (int)strlen(str));
  //printBacktrace();
  fflush(stdout);
  return result;
}



ssize_t __wrap_write(int fd, const void *buf, size_t count) {
  int prev_inside_wrapper = inside_wrapper;
  inside_wrapper = 1;

  ssize_t result = __real_write(fd, buf, count);
  
  if (!prev_inside_wrapper) {
    printf("  call to write, fd=%d, count=%d\n", fd, (int)count);
    printBacktrace();
    /* write(1, buf, count); */
    fflush(stdout);
  }

  inside_wrapper = prev_inside_wrapper;

  return result;
}


ssize_t __wrap_pwrite(int fd, const void *buf, size_t count, off_t offset) {
  int prev_inside_wrapper = inside_wrapper;
  inside_wrapper = 1;

  ssize_t result = __real_pwrite(fd, buf, count, offset);
  
  if (!prev_inside_wrapper) {
    printf("  call to pwrite, fd=%d, count=%d, offset=%ld\n",
             fd, (int)count, (long)offset);
    printBacktrace();
    /* write(1, buf, count); */
    fflush(stdout);
  }

  inside_wrapper = prev_inside_wrapper;

  return result;
}


ssize_t __wrap_writev(int fd, const struct iovec *iov, int iovcnt) {
  int prev_inside_wrapper = inside_wrapper;
  inside_wrapper = 1;

  ssize_t result = __real_writev(fd, iov, iovcnt);
  
  if (!prev_inside_wrapper) {
    printf("  call to writev, fd=%d, iovcnt=%d\n", fd, iovcnt);
    printBacktrace();
    fflush(stdout);
  }

  inside_wrapper = prev_inside_wrapper;

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

