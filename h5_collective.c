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
#define MIMIC_HDF5_DATATYPE 1
#define SET_MPI_ATOMIC 0

int rank, np;

/* note: these are incompatible with Darshan. If these are enabled, run
   "module unload darshan" before compiling. */
#include "wrapper_fns.c"

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
    /* printf("MPI_ATOMIC=%d\n", SET_MPI_ATOMIC); */
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
      p.data[y * p.col_count + x] = rank + .01*y + .0001 * (x + p.col_start);
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
      continue;
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
      continue;
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


/*
  for 10x10 data, here is the filetype (in the call to MPI_File_set_view
    vector count=1, blocklength=10, stride=1 of:
      create_resized lb=0, extent=80 of:
        create_hindexed count=1 (1 at offset 0) of:
          vector count=1, blocklength=5, stride=1 of:
            contiguous 8 MPI_BYTE
*/
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

  H5Fset_mpi_atomicity(file_id, SET_MPI_ATOMIC);

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

#if MIMIC_HDF5_DATATYPE
  /* define the full array on disk */
  MPI_Datatype double_bytes;
  MPI_Type_contiguous(sizeof(double), MPI_BYTE, &double_bytes);
  MPI_Type_commit(&double_bytes);
  element_type = double_bytes;
  MPI_Datatype type1_vec, type2_hindex, type3_resized;
  MPI_Type_vector(1, p->col_count, 1, element_type, &type1_vec);
  int hindex_len = 1;
  MPI_Aint hindex_offset = p->col_start * sizeof(double);
  MPI_Type_create_hindexed(1, &hindex_len, &hindex_offset, type1_vec, &type2_hindex);
  MPI_Type_free(&type1_vec);
  MPI_Type_create_resized(type2_hindex, 0, sizeof(double) * p->cols, &type3_resized);
  MPI_Type_free(&type2_hindex);
  MPI_Type_vector(1, p->rows, 1, type3_resized, &file_type);
  MPI_Type_free(&type3_resized);

  /* define the subarray in memory */
  MPI_Type_vector(1, p->rows * p->col_count, 1, element_type, &memory_type);

#else
  /* define the full array on disk */
  int full_sizes[2] = {p->rows, p->cols};
  int part_sizes[2] = {p->rows, p->col_count};
  int part_starts[2] = {0, p->col_start};
  MPI_Type_create_subarray(2, full_sizes, part_sizes, part_starts,
                           MPI_ORDER_C, element_type, &file_type);

  /* define the subarray in memory */
  MPI_Type_contiguous(p->rows * p->col_count, element_type, &memory_type);

#endif

  MPI_Type_commit(&file_type);
  MPI_Type_commit(&memory_type);

  /*
  if (rank == 0) {
    char *s = mpiTypeSignature(file_type, NULL);
    printf("mpiTypeSignature(file_type) = \"%s\"\n", s);
    free(s);
  }
  */


  /* open the file */
  remove(filename);
  status = MPI_File_open(comm, filename, MPI_MODE_RDWR | MPI_MODE_CREATE,
                         p->info, &file);
  if (status != MPI_SUCCESS) {
    printf("Failed to open \"%s\". Cannot do MPI-IO test.\n", filename);
    return;
  }

  MPI_File_set_atomicity(file, SET_MPI_ATOMIC);
  
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
