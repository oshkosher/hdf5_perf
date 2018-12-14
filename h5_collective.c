#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <hdf5.h>

/*
  Profile MPI_File_open() and MPI_File_set_view() calls.
  Check the info values before and after the call.
  Check which form of 'write' HDF5 is using.
*/

#define DATASETNAME "/mydata"
#define PRINT_PARTITIONS 0

#define PROFILE_FNS 0

int rank, np;

typedef struct {
  int rows, cols;  /* global size */
  int col_start, col_count;  /* this process' data */
  MPI_Info info;  /* striping settings */
  double *data;  /* rows*col_count entries */
} Params;

int parseArgs(int argc, char **argv, const char **filename,
              int *rows, int *cols,
              int *stripe_count, int *stripe_len);
void makeFilenames(const char *base_filename, char **filename_h5,
                   char **fliename_mpiio);
void printHelp();
void partition(int total_size, int np, int rank, 
               int *offset, int *size);
int check_internal(int status, char *filename, int line_no);
void writeHDF5File(const char *filename, Params *p);
void writeMPIIOFile(const char *filename, Params *p);

/* profiling */
int MPI_File_open
  (MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File *fh);
int PMPI_File_open
  (MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File *fh);
int MPI_File_set_view
  (MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,
   const char *datarep, MPI_Info info);
int PMPI_File_set_view
  (MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,
   const char *datarep, MPI_Info info);

#define check(s) check_internal(s, __FILE__, __LINE__)

int main(int argc, char **argv) {
  double timer;
  int stripe_count = -1, stripe_len = -1;
  const char *base_filename = NULL;
  char *filename_h5 = NULL, *filename_mpiio = NULL;

  Params p;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (parseArgs(argc, argv, &base_filename, &p.rows, &p.cols,
                &stripe_count, &stripe_len))
    exit(1);


  if (rank==0) {
    double mb = sizeof(double) * p.rows * p.cols / (1024.0 * 1024);
    printf("Write %dx%d grid (%.1f MiB) with %d ranks.\n",
           p.rows, p.cols, mb, np);
  }

  makeFilenames(base_filename, &filename_h5, &filename_mpiio);
  
  /* each process writes every row, but only a portion of each row */
  partition(p.cols, np, rank, &p.col_start, &p.col_count);

#if PRINT_PARTITIONS
  for (int r=0; r < np; r++) {
    if (rank == r)
      printf("[%d] %d..%d\n", rank, p.col_start,
             p.col_start + p.col_count - 1);
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  timer = MPI_Wtime();
  size_t data_size_bytes = sizeof(double) * p.rows * p.col_count;
  p.data = (double*) malloc(data_size_bytes);
  if (!p.data) {
    printf("[%d] ERROR failed to allocate %ld bytes\n",
           rank, (long)(data_size_bytes));
    exit(1);
  }
  for (int y=0; y < p.rows; y++) {
    for (int x=0; x < p.col_count; x++) {
      p.data[y * p.col_count + x] = rank + .01*y + .0001 * x;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0)
    printf("Init time: %.6f\n", MPI_Wtime() - timer);

  /* set the file striping parameters */
  MPI_Info_create(&p.info);
  char tmp_str[50];
  sprintf(tmp_str, "%d", stripe_count);
  MPI_Info_set(p.info, "striping_factor", tmp_str);
  sprintf(tmp_str, "%d", stripe_len);
  MPI_Info_set(p.info, "striping_unit", tmp_str);

  writeHDF5File(filename_h5, &p);
  writeMPIIOFile(filename_mpiio, &p);

  MPI_Info_free(&p.info);
  free(filename_h5);
  free(filename_mpiio);
  free(p.data);

  MPI_Finalize();

  return 0;
}


int parseArgs(int argc, char **argv, const char **filename,
              int *rows, int *cols,
              int *stripe_count, int *stripe_len) {

  int argno = 1;
  *stripe_count = 4;
  *stripe_len = 1 * (1 << 20);

  while (argno < argc && argv[argno][0] == '-') {
    const char *arg = argv[argno];
    
    if (!strcmp(arg, "-stripe_count")) {
      argno++;
      if (argno >= argc) printHelp();
      arg = argv[argno];
      if (1 != sscanf(arg, "%d", stripe_count) ||
          *stripe_count < 0) {
        printf("Invalid stripe count: \"%s\"\n", arg);
        return -1;
      }
      argno++;
    }
    
    if (!strcmp(arg, "-stripe_len")) {
      argno++;
      if (argno >= argc) printHelp();
      arg = argv[argno];
      if (1 != sscanf(arg, "%d", stripe_len) ||
          *stripe_len < 0) {
        printf("Invalid stripe len: \"%s\"\n", arg);
        return -1;
      }
      *stripe_len *= (1<<20);
      argno++;
    }

    else
      printHelp();
  }

  if (argc - argno != 3) printHelp();

  *filename = argv[argno++];
  
  if (1 != sscanf(argv[argno], "%d", rows) || *rows <= 0) {
    fprintf(stderr, "Invalid #rows: %s\n", argv[argno]);
    return -1;
  }
  argno++;
  
  if (1 != sscanf(argv[argno], "%d", cols) || *cols <= 0) {
    fprintf(stderr, "Invalid #cols: %s\n", argv[argno]);
    return -1;
  }
  argno++;

  return 0;
}


void makeFilenames(const char *base_filename, char **filename_h5,
                   char **filename_mpiio) {
  *filename_h5 = (char*) malloc(strlen(base_filename) + 4);
  sprintf(*filename_h5, "%s.h5", base_filename);

  *filename_mpiio = (char*) malloc(strlen(base_filename) + 7);
  sprintf(*filename_mpiio, "%s.mpiio", base_filename);
}


void printHelp() {
  fprintf(stderr, "\n"
          "  h5_collective [opt] <filename> <#rows> <#cols>\n"
          "  opt:\n"
          "    -stripe_count <n> : number of file stripes\n"
          "    -stripe_len <n> : length of file stripes, in MiB\n"
          "  Compare performance of HDF5 collective IO vs. MPI-IO direct.\n"
          "\n");
  exit(1);
}


void partition(int total_size, int np, int rank, 
               int *offset, int *size) {
  *offset = rank * total_size / np;
  int end = (rank + 1) * total_size / np;
  *size = end - *offset;
}


int check_internal(int status, char *filename, int line_no) {
  if (status < 0)
    printf("[%d] %s:%d status %d\n", rank, filename, line_no, status);
  return status;
}


void writeHDF5File(const char *filename, Params *p) {
  hid_t file_id=-1, dataset_id=-1, dataspace_id=-1, memspace_id=-1;
  int status=-1;
  double timer_open, timer_write, timer_close;
  MPI_Comm comm = MPI_COMM_WORLD;

  hid_t file_create_properties = H5P_DEFAULT;
  hid_t file_access_properties = H5P_DEFAULT;

  /* open the file */
  MPI_Barrier(comm);
  timer_open = MPI_Wtime();

  file_access_properties = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(file_access_properties, comm, p->info);

  remove(filename);

  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, 
                      file_create_properties, file_access_properties);
  if (file_id < 0) {
    fprintf(stderr, "Failed to create %s.\n", filename);
    exit(1);
  }

  H5Pclose(file_access_properties);

  /* create dataspace */
  hsize_t full_dims[2] = {p->rows, p->cols};
  dataspace_id = H5Screate_simple(2, full_dims, NULL);
  check(dataspace_id);

  /* create the dataset */
  dataset_id = H5Dcreate2
    (file_id, DATASETNAME, H5T_NATIVE_DOUBLE, dataspace_id,
     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  check(dataset_id);

  /* define the subset of data I'll write */
  hsize_t file_offsets[2] = {0, p->col_start};
  hsize_t file_sizes[2] = {p->rows, p->col_count};
  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, file_offsets,
                               NULL, file_sizes, NULL);
  check(status);
  
  /* define the shape of the memory buffer */
  hsize_t mem_offsets[2] = {0, 0};
  hsize_t mem_sizes[2] = {p->rows, p->col_count};
  memspace_id = H5Screate_simple(2, mem_sizes, NULL);
  check(memspace_id);

  status = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mem_offsets,
                               NULL, mem_sizes, NULL);
  check(status);

  MPI_Barrier(MPI_COMM_WORLD);
  timer_open = MPI_Wtime() - timer_open;

  /* write content */
  timer_write = MPI_Wtime();
  hid_t io_prop = H5Pcreate(H5P_DATASET_XFER);
  status = H5Pset_dxpl_mpio(io_prop, H5FD_MPIO_COLLECTIVE);
  check(status);

  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,
                    io_prop, p->data);
  check(status);
  H5Pclose(io_prop);
  MPI_Barrier(MPI_COMM_WORLD);
  timer_write = MPI_Wtime() - timer_write;

  /* close the file */  
  timer_close = MPI_Wtime();
  check(H5Sclose(memspace_id));
  check(H5Dclose(dataset_id));
  check(H5Sclose(dataspace_id));
  check(H5Fclose(file_id));
  MPI_Barrier(MPI_COMM_WORLD);
  timer_close = MPI_Wtime() - timer_close;

  if (rank == 0)
    printf("HDF5 write %s\n"
           "  Open %.3fs, write %.3fs, close %.3fs. Write %.1f MiB/s\n",
           filename,
           timer_open, timer_write, timer_close,
           (p->rows*p->cols*sizeof(double) / (timer_write * (1<<20))));
}


void writeMPIIOFile(const char *filename, Params *p) {
  int rank, status;
  MPI_File file;
  MPI_Datatype file_type, memory_type;
  MPI_Comm comm = MPI_COMM_WORLD;
  double timer_open, timer_write, timer_close;

  MPI_Comm_rank(comm, &rank);

  timer_open = MPI_Wtime();

  MPI_Datatype element_type = MPI_DOUBLE;

  MPI_Datatype double_bytes;
  MPI_Type_contiguous(sizeof(double), MPI_BYTE, &double_bytes);
  MPI_Type_commit(&double_bytes);

  element_type = double_bytes;
  // element_type = MPI_DOUBLE;

  /* define the full array on disk */
  int full_sizes[2] = {p->rows, p->cols};
  int part_sizes[2] = {p->rows, p->col_count};
  int part_starts[2] = {0, p->col_start};
  MPI_Type_create_subarray(2, full_sizes, part_sizes, part_starts,
                           MPI_ORDER_C, element_type, &file_type);
  MPI_Type_commit(&file_type);

  /* define the subarray in memory */
  MPI_Type_contiguous(p->rows * p->col_count, element_type, &memory_type);
  MPI_Type_commit(&memory_type);

  /* open the file */
  remove(filename);
  status = MPI_File_open(comm, filename, MPI_MODE_RDWR | MPI_MODE_CREATE,
                         p->info, &file);
  if (status != MPI_SUCCESS) {
    printf("Failed to open \"%s\". Cannot do MPI-IO test.\n", filename);
    return;
  }
  
  /* set the file view */
  MPI_File_set_view(file, 0, MPI_BYTE, file_type, "native", p->info);

  MPI_Barrier(comm);
  timer_open = MPI_Wtime() - timer_open;

  /* write the data */
  MPI_Status status2;
  int count_written;
  timer_write = MPI_Wtime();
  status = MPI_File_write_at_all(file, 0, p->data, 1, memory_type, &status2);
  if (status != MPI_SUCCESS) {
    printf("[%d] MPI_File_write failed, error = %d\n", rank, status);
    return;
  }
  MPI_Get_count(&status2, memory_type, &count_written);
  if (count_written != 1)
    printf("[%d] MPI_File_write_at_all count=%d\n", rank, count_written);
  MPI_Barrier(comm);
  timer_write = MPI_Wtime() - timer_write;

  timer_close = MPI_Wtime();
  if (element_type != MPI_DOUBLE)
    MPI_Type_free(&element_type);
  MPI_Type_free(&memory_type);
  MPI_Type_free(&file_type);
  MPI_File_close(&file);
  timer_close = MPI_Wtime() - timer_close;

  MPI_Barrier(comm);
  if (rank == 0)
    printf("MPI-IO write %s\n"
           "  Open %.3fs, write %.3fs, close %.3fs. Write %.1f MiB/s\n",
           filename,
           timer_open, timer_write, timer_close,
           (p->rows*p->cols*sizeof(double) / (timer_write * (1<<20))));
}


static int amode_flags[] = {
  MPI_MODE_CREATE,
  MPI_MODE_RDONLY,
  MPI_MODE_WRONLY,
  MPI_MODE_RDWR,
  MPI_MODE_DELETE_ON_CLOSE,
  MPI_MODE_UNIQUE_OPEN,
  MPI_MODE_EXCL,
  MPI_MODE_APPEND,
  MPI_MODE_SEQUENTIAL
};

static char* amode_names[] = {
  "MPI_MODE_CREATE",
  "MPI_MODE_RDONLY",
  "MPI_MODE_WRONLY",
  "MPI_MODE_RDWR",
  "MPI_MODE_DELETE_ON_CLOSE",
  "MPI_MODE_UNIQUE_OPEN",
  "MPI_MODE_EXCL",
  "MPI_MODE_APPEND",
  "MPI_MODE_SEQUENTIAL"
};  

static char *formatAmode(char buf[174], int amode) {
  /* concatenation of all 9 possible modes with 3-byte separators
     = 173 chars */
  assert(sizeof(amode_flags) / sizeof(amode_flags[0]) == 9);
  char *p = buf;
  *p = 0;
  for (int i=0; i < 9; i++) {
    if (amode & amode_flags[i]) {
      if (p != buf) {
        strcat(p, " | ");
        p += 3;
      }
      strcat(p, amode_names[i]);
      p += strlen(amode_names[i]);
    }
  }
  return buf;
}


static void printInfo(MPI_Info info, const char *prefix) {
  char key[MPI_MAX_INFO_KEY + 1], *value;
  int value_buf_len = 256;

  if (info == MPI_INFO_NULL) {
    printf("%s  MPI_INFO_NULL\n", prefix);
    return;
  }

  value = (char*) malloc(value_buf_len);
  assert(value);

  int nkeys;
  MPI_Info_get_nkeys(info, &nkeys);
  for (int keyno=0; keyno < nkeys; keyno++) {
    MPI_Info_get_nthkey(info, keyno, key);
    int value_len, key_defined;
    MPI_Info_get_valuelen(info, key, &value_len, &key_defined);
    assert(key_defined);

    // make sure the buffer is large enough for the value
    if (value_len+1 > value_buf_len) {
      value_buf_len *= 2;
      if (value_len+1 > value_buf_len)
        value_buf_len = value_len+1;
      value = (char*) realloc(value, value_buf_len);
      assert(value);
    }

    MPI_Info_get(info, key, value_buf_len-1, value, &key_defined);
    assert(key_defined);

    printf("%s  MPI_Info \"%s\": \"%s\"\n", prefix, key, value);
  }
  free(value);
}


void printFileInfo(MPI_File f, const char *prefix) {
  MPI_Info info;

  MPI_File_get_info(f, &info);

  printInfo(info, prefix);
  
  MPI_Info_free(&info);
}


#if PROFILE_FNS

/* return a string representing this datatype. If it's a primitive type, return
   the name (MPI_DOUBLE, MPI_INT). Otherwise return a hex string. */
static char* mpiTypeTostr(char outstr[22], MPI_Datatype dt) {
  switch (dt) {
  case MPI_CHAR: strcpy(outstr, "MPI_CHAR"); break;
  case MPI_SIGNED_CHAR: strcpy(outstr, "MPI_SIGNED_CHAR"); break;
  case MPI_UNSIGNED_CHAR: strcpy(outstr, "MPI_UNSIGNED_CHAR"); break;
  case MPI_BYTE: strcpy(outstr, "MPI_BYTE"); break;
  case MPI_WCHAR: strcpy(outstr, "MPI_WCHAR"); break;
  case MPI_SHORT: strcpy(outstr, "MPI_SHORT"); break;
  case MPI_UNSIGNED_SHORT: strcpy(outstr, "MPI_UNSIGNED_SHORT"); break;
  case MPI_INT: strcpy(outstr, "MPI_INT"); break;
  case MPI_UNSIGNED: strcpy(outstr, "MPI_UNSIGNED"); break;
  case MPI_LONG: strcpy(outstr, "MPI_LONG"); break;
  case MPI_UNSIGNED_LONG: strcpy(outstr, "MPI_UNSIGNED_LONG"); break;
  case MPI_FLOAT: strcpy(outstr, "MPI_FLOAT"); break;
  case MPI_DOUBLE: strcpy(outstr, "MPI_DOUBLE"); break;
  case MPI_LONG_DOUBLE: strcpy(outstr, "MPI_LONG_DOUBLE"); break;
  case MPI_LONG_LONG_INT: strcpy(outstr, "MPI_LONG_LONG_INT"); break;
  case MPI_UNSIGNED_LONG_LONG: strcpy(outstr, "MPI_UNSIGNED_LONG_LONG"); break;
  default:
    sprintf(outstr, "0x%08x", (int)dt);
  }
  return outstr;
}


int MPI_File_open
  (MPI_Comm comm, const char *filename, int amode, MPI_Info info,
   MPI_File *fh) {
  const char *prefix = "[0]  ";

  /*
  if (rank==0) {
    printf("[%d] MPI_File_open %s info:\n", rank, filename);
    printInfo(info, prefix);
  }
  */

  int result = PMPI_File_open(comm, filename, amode, info, fh);
  
  if (rank==0) {
    printf("[%d] MPI_File_open(\"%s\", fh=%p)\n", rank, filename, (void*)*fh);
    /* printFileInfo(*fh, prefix); */
  }
  
  return result;
}


int MPI_File_set_view
  (MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,
   const char *datarep, MPI_Info info) {
  const char *prefix = "[0]  ";
  
  if (rank==0) {
    char buf1[50], buf2[50];
    printf("[%d] MPI_File_set_view fh=%p disp=%ld, etype=%s, filetype=%s, datarep=%s\n",
           rank, (void*)fh, (long)disp, mpiTypeTostr(buf1, etype), mpiTypeTostr(buf2, filetype), datarep);
    /* printInfo(info, prefix); */
  }

  return PMPI_File_set_view(fh, disp, etype, filetype, datarep, info);
}



int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void * buf, int count, MPI_Datatype datatype, MPI_Status * status) {
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_File_write_at(fh=%p, offset=%ld, count=%d, type=%s)\n",
           (void*)fh, (long)offset, count, mpiTypeTostr(buf, datatype));
  }
  return PMPI_File_write_at(fh, offset, buf, count, datatype, status);
}

int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void * buf, int count, MPI_Datatype datatype, MPI_Status * status) {
  if (rank==0) printf("[0] MPI_File_write_at_all(file=%p, offset=%ld, count=%d, type=0x%08x)\n", (void*)fh, (long)offset, count, (int)datatype);
  return PMPI_File_write_at_all(fh, offset, buf, count, datatype, status);
}

int MPI_File_iwrite_at(MPI_File fh, MPI_Offset offset, const void * buf, int count, MPI_Datatype datatype, MPIO_Request * request) {
  if (rank==0) printf("[0] MPI_File_iwrite_at(%p)\n", (void*)fh);
  return PMPI_File_iwrite_at(fh, offset, buf, count, datatype, request);
}

int MPI_File_write(MPI_File fh, const void * buf, int count, MPI_Datatype datatype, MPI_Status * status) {
  if (rank==0) printf("[0] MPI_File_write(%p)\n", (void*)fh);
  return PMPI_File_write(fh, buf, count, datatype, status);
}

int MPI_File_write_all(MPI_File fh, const void * buf, int count, MPI_Datatype datatype, MPI_Status * status) {
  if (rank==0) printf("[0] MPI_File_write_all(%p)\n", (void*)fh);
  return PMPI_File_write_all(fh, buf, count, datatype, status);
}

int MPI_File_iwrite(MPI_File fh, const void * buf, int count, MPI_Datatype datatype, MPIO_Request * request) {
  if (rank==0) printf("[0] MPI_File_iwrite(%p)\n", (void*)fh);
  return PMPI_File_iwrite(fh, buf, count, datatype, request);
}

int MPI_File_write_shared(MPI_File fh, const void * buf, int count, MPI_Datatype datatype, MPI_Status * status) {
  if (rank==0) printf("[0] MPI_File_write_shared(%p)\n", (void*)fh);
  return PMPI_File_write_shared(fh, buf, count, datatype, status);
}

int MPI_File_iwrite_shared(MPI_File fh, const void * buf, int count, MPI_Datatype datatype, MPIO_Request * request) {
  if (rank==0) printf("[0] MPI_File_iwrite_shared(%p)\n", (void*)fh);
  return PMPI_File_iwrite_shared(fh, buf, count, datatype, request);
}

int MPI_File_write_ordered(MPI_File fh, const void * buf, int count, MPI_Datatype datatype, MPI_Status * status) {
  if (rank==0) printf("[0] MPI_File_write_ordered(%p)\n", (void*)fh);
  return PMPI_File_write_ordered(fh, buf, count, datatype, status);
}

int MPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, const void * buf, int count, MPI_Datatype datatype) {
  if (rank==0) printf("[0] MPI_File_write_at_all_begin(%p)\n", (void*)fh);
  return PMPI_File_write_at_all_begin(fh, offset, buf, count, datatype);
}

int MPI_File_write_at_all_end(MPI_File fh, const void * buf, MPI_Status * status) {
  if (rank==0) printf("[0] MPI_File_write_at_all_end(%p)\n", (void*)fh);
  return PMPI_File_write_at_all_end(fh, buf, status);
}

int MPI_File_write_all_begin(MPI_File fh, const void * buf, int count, MPI_Datatype datatype) {
  if (rank==0) printf("[0] MPI_File_write_all_begin(%p)\n", (void*)fh);
  return PMPI_File_write_all_begin(fh, buf, count, datatype);
}

int MPI_File_write_all_end(MPI_File fh, const void * buf, MPI_Status * status) {
  if (rank==0) printf("[0] MPI_File_write_all_end(%p)\n", (void*)fh);
  return PMPI_File_write_all_end(fh, buf, status);
}

int MPI_File_write_ordered_begin(MPI_File fh, const void * buf, int count, MPI_Datatype datatype) {
  if (rank==0) printf("[0] MPI_File_write_ordered_begin(%p)\n", (void*)fh);
  return PMPI_File_write_ordered_begin(fh, buf, count, datatype);
}

int MPI_File_write_ordered_end(MPI_File fh, const void * buf, MPI_Status * status) {
  if (rank==0) printf("[0] MPI_File_write_ordered_end(%p)\n", (void*)fh);
  return PMPI_File_write_ordered_end(fh, buf, status);
}

int MPI_File_iwrite_at_all(MPI_File fh, MPI_Offset offset, const void * buf, int count, MPI_Datatype datatype, MPI_Request * request) {
  if (rank==0) printf("[0] MPI_File_iwrite_at_all(%p)\n", (void*)fh);
  return PMPI_File_iwrite_at_all(fh, offset, buf, count, datatype, request);
}

int MPI_File_iwrite_all(MPI_File fh, const void * buf, int count, MPI_Datatype datatype, MPI_Request * request) {
  if (rank==0) printf("[0] MPI_File_iwrite_all(%p)\n", (void*)fh);
  return PMPI_File_iwrite_all(fh, buf, count, datatype, request);
}


int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype * newtype) {
  int result = PMPI_Type_contiguous(count, oldtype, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_contiguous(count=%d, oldtype=%s, newtype=0x%08x)\n",
           count, mpiTypeTostr(buf, oldtype), *(int*)newtype);
  }
  return result;
}

int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype * newtype) {
  int result = PMPI_Type_vector(count, blocklength, stride, oldtype, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_vector(count=%d, blocklength=%d, stride=%d, oldtype=%s, newtype=0x%08x)\n",
           count, blocklength, stride, mpiTypeTostr(buf, oldtype), *(int*)newtype);
  }
  return result;
}

int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype * newtype) {
  int result = PMPI_Type_hvector(count, blocklength, stride, oldtype, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_hvector(%d, %d, %ld, %s, 0x%08x)\n",
           count, blocklength, stride, mpiTypeTostr(buf, oldtype), *(int*)newtype);
  }
  return result;
}

int MPI_Type_indexed(int count, const int * array_of_blocklengths, const int * array_of_displacements, MPI_Datatype oldtype, MPI_Datatype * newtype) {
  int result = PMPI_Type_indexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_indexed(%d, ..., ..., %s, 0x%08x)\n",
           count, mpiTypeTostr(buf, oldtype), *(int*)newtype);
  }
  return result;
}

int MPI_Type_hindexed(int count, int * array_of_blocklengths, MPI_Aint * array_of_displacements, MPI_Datatype oldtype, MPI_Datatype * newtype) {
  int result = PMPI_Type_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_hindexed(%d, ..., ..., %s, 0x%08x)\n",
           count, mpiTypeTostr(buf, oldtype), *(int*)newtype);
  }
  return result;
}

int MPI_Type_struct(int count, int * array_of_blocklengths, MPI_Aint * array_of_displacements, MPI_Datatype * array_of_types, MPI_Datatype * newtype) {
  int result = PMPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_struct(%d, ..., ..., 0x%08x)\n",
           count, *(int*)newtype);
  }
  return result;
}


int MPI_Type_create_darray(int size , int rankarg , int ndims , const int array_of_gsizes [], const int array_of_distribs [], const int array_of_dargs [], const int array_of_psizes [], int order , MPI_Datatype oldtype , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_darray(size, rankarg, ndims, array_of_gsizes, array_of_distribs, array_of_dargs, array_of_psizes, order, oldtype, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_darray(... %s, 0x%08x)\n",
           mpiTypeTostr(buf, oldtype), *(int*)newtype);
  }
  return result;
}

int MPI_Type_create_hindexed(int count , const int array_of_blocklengths [], const MPI_Aint array_of_displacements [], MPI_Datatype oldtype , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_hindexed(count=%d, oldtype=%s, newtype=0x%08x)\n",
           count, mpiTypeTostr(buf, oldtype), *(int*)newtype);
    /*
    for (int i=0; i < count; i++) {
      if (i > 0) printf(", ");
      printf("%d@%ld", array_of_blocklengths[i], (long)array_of_displacements[i]);
    }
    putchar('\n');
    */
  }
  return result;
}

int MPI_Type_create_hvector(int count , int blocklength , MPI_Aint stride , MPI_Datatype oldtype , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_hvector(count, blocklength, stride, oldtype, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_hvector(... %s, 0x%08x)\n",
           mpiTypeTostr(buf, oldtype), *(int*)newtype);
  }
  return result;
}

int MPI_Type_create_indexed_block(int count , int blocklength , const int array_of_displacements [], MPI_Datatype oldtype , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_indexed_block(count, blocklength, array_of_displacements, oldtype, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_indexed_block(... %s, 0x%08x)\n",
           mpiTypeTostr(buf, oldtype), *(int*)newtype);
  }
  return result;
}

int MPI_Type_create_hindexed_block(int count , int blocklength , const MPI_Aint array_of_displacements [], MPI_Datatype oldtype , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_hindexed_block(count, blocklength, array_of_displacements, oldtype, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_hindexed_block(... %s, 0x%08x)\n",
           mpiTypeTostr(buf, oldtype), *(int*)newtype);
  }
  return result;
}

int MPI_Type_create_resized(MPI_Datatype oldtype , MPI_Aint lb , MPI_Aint extent , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_resized(oldtype, lb, extent, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_resized(lb=%ld, extent=%ld, oldtype=%s, newtype=0x%08x)\n",
           (long)lb, (long)extent, mpiTypeTostr(buf, oldtype), *(int*)newtype);
  }
  return result;
}

int MPI_Type_create_struct(int count , const int array_of_blocklengths [], const MPI_Aint array_of_displacements [], const MPI_Datatype array_of_types [], MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
  if (rank==0) {
    char buf[50];
    /*
    printf("[0] MPI_Type_create_struct(... newtype=0x%08x)\n",
           *(int*)newtype);
    */

    printf("[0] MPI_Type_create_struct(count=%d, newtype=0x%08x) {", count, *(int*)newtype);
    for (int i=0; i < count; i++) {
      char buf[50];
      if (i > 0) printf(", ");
      printf("%d %s @ %ld", array_of_blocklengths[i], mpiTypeTostr(buf, array_of_types[i]), (long)array_of_displacements[i]);
    }
    printf("}\n");
  }
  return result;
}

int MPI_Type_create_subarray(int ndims , const int array_of_sizes [], const int array_of_subsizes [], const int array_of_starts [], int order , MPI_Datatype oldtype , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, order, oldtype, newtype);
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_subarray(... %s, 0x%08x)\n",
           mpiTypeTostr(buf, oldtype), *(int*)newtype);
  }
  return result;
}

int MPI_Type_commit(MPI_Datatype *datatype) {
  int result = PMPI_Type_commit(datatype);
  if (rank==0) {
    printf("[0] MPI_Type_commit(0x%08x)\n", *datatype);
  }
  return result;
}


int MPI_Type_free(MPI_Datatype *datatype) {
  if (rank==0) {
    printf("[0] MPI_Type_free(0x%08x)\n", *datatype);
  }
  return PMPI_Type_free(datatype);
}


#endif
