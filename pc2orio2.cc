#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <string>
#include <vector>
#include <mpi.h>
#include <hdf5.h>

using std::string;
using std::vector;

#ifndef H5_HAVE_PARALLEL
#error No HDF5 parallel
#endif

#define PRINT_TIMESTAMPS 0

/*
  Test the writing of a file that contains multiple 3-d grids split equally
  across ranks.

  Compare HDF5 where each grid is a dataspace and MPI-IO, where each grid
  is appended to the file.

  Ed Karrels, edk@illinois.edu, December 2018
*/


int rank, np;   // rank & size of MPI_COMM_WORLD
double time0;


class Triple {
  // save the result of str(), so it doesn't have to be stored in a temp
  // variable to pass it to printf().
  mutable char str_buf[50];

public:
  int coord[3];

  int& operator[](size_t i) {return coord[i];}
  int operator[](size_t i) const {return coord[i];}

  const int* array() const {return &coord[0];}

  // parse a string in the form "x,y,z". Return false on error.
  bool parse(const char *str) {
    if (3 == sscanf(str, "%d,%d,%d", coord, coord+1, coord+2))
      return true;
    else
      return false;
  }

  int product() const {
    return coord[0] * coord[1] * coord[2];
  }

  // format as a string in the form "x,y,z"
  const char *str() const {
    sprintf(str_buf, "%d,%d,%d", coord[0], coord[1], coord[2]);
    return str_buf;
  }

  // return true iff all 3 coordinates are positive
  bool isPositive() const {
    return coord[0] > 0 && coord[1] > 0 && coord[2] > 0;
  }

  void toHsize(hsize_t dest[3]) {
    for (int i=0; i < 3; i++) dest[i] = coord[i];
  }
};


struct Params {
  Triple global_size;  /* global size */
  Triple local_starts, local_counts;  /* this process's subgrid */
  int grid_count;  /* how many grids to write */

  MPI_Comm cart_comm;
  int cart_rank;  // rank in cart_comm
  Triple cart_coords;  // 3-d coordinates

  vector<double> data;  /* this process's data */
  MPI_Info info;  /* striping settings */
};


bool parseArgs(int argc, char **argv, string &filename,
               Triple &global_size, Triple &rank_splits,
               int &grid_count,
               int &stripe_count, int &stripe_len);
void fail();
void printHelp();
bool setupData(Params &p, const Triple &rank_splits,
               int stripe_count, int stripe_len);
void partition(int total_size, int np, int rank, 
               int &offset, int &size);
void timestamp(const char *name);
int h5_check_internal(int status, const char *filename, int line_no);

void writeHDF5File(const char *filename, Params &p);
void writeMPIIOFile(const char *filename, Params &p,
                    bool use_vector_type = false);

#define h5check(s) h5_check_internal(s, __FILE__, __LINE__)


int main(int argc, char **argv) {
  double timer;
  int stripe_count = -1, stripe_len = -1;
  string base_filename, filename_h5, filename_mpiio;
  Triple rank_splits;
  Params p;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  time0 = MPI_Wtime();
  timestamp("after MPI_Init");

  // printf("Hello from rank %d\n", rank);

  if (!parseArgs(argc, argv, base_filename, p.global_size,
                 rank_splits, p.grid_count, stripe_count, stripe_len))
    fail();

  if (rank==0) {
    double mb = sizeof(double) * p.global_size.product() * p.grid_count
      / (1024.0 * 1024);
    printf("Write %dx(%s) grids (%.1f MiB), %d ranks (%s), "
           "stripe_count=%d, stripe_len=%d\n",
           p.grid_count, p.global_size.str(),
           mb, np, rank_splits.str(),
           stripe_count, stripe_len);
  }

  filename_h5 = base_filename + ".h5";
  filename_mpiio = base_filename + ".mpiio";

  MPI_Barrier(MPI_COMM_WORLD);
  timer = MPI_Wtime();

  timestamp("before setupData");
  if (!setupData(p, rank_splits, stripe_count, stripe_len))
    fail();
  timestamp("after setupData");

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0)
    printf("Init time: %.6f\n", MPI_Wtime() - timer);

  timestamp("before writeHDF5File");
  writeHDF5File(filename_h5.c_str(), p);
  timestamp("after writeHDF5File");
  writeMPIIOFile(filename_mpiio.c_str(), p);
  timestamp("after writeMPIIOFile");

  MPI_Info_free(&p.info);
  MPI_Comm_free(&p.cart_comm);
  MPI_Finalize();

  return 0;
}


bool parseArgs(int argc, char **argv, string &filename,
               Triple &global_size, Triple &rank_splits,
               int &grid_count,
               int &stripe_count, int &stripe_len) {

  int argno = 1;

  // default 4x1MB stripes
  stripe_count = 4;
  stripe_len = 1 * (1 << 20);

  while (argno < argc && argv[argno][0] == '-') {
    const char *arg = argv[argno];
    
    if (!strcmp(arg, "-stripe_count")) {
      argno++;
      if (argno >= argc) printHelp();
      arg = argv[argno];
      if (1 != sscanf(arg, "%d", &stripe_count) ||
          stripe_count < 0) {
        if (rank==0)
          printf("Invalid stripe count: \"%s\"\n", arg);
        return false;
      }
      argno++;
      continue;
    }
    
    if (!strcmp(arg, "-stripe_len")) {
      argno++;
      if (argno >= argc) printHelp();
      arg = argv[argno];
      if (1 != sscanf(arg, "%d", &stripe_len) ||
          stripe_len < 0) {
        if (rank==0)
          printf("Invalid stripe len: \"%s\"\n", arg);
        return false;
      }
      stripe_len *= (1<<20);
      argno++;
      continue;
    }

    else
      printHelp();
  }

  if (argc - argno != 4) printHelp();

  filename = argv[argno++];

  if (!global_size.parse(argv[argno]) ||
      !global_size.isPositive()) {
    if (rank == 0)
      printf("Invalid global size: %s\n", argv[argno]);
    return false;
  }
  argno++;

  if (!rank_splits.parse(argv[argno]) ||
      !rank_splits.isPositive()) {
    if (rank == 0)
      printf("Invalid global size: %s\n", argv[argno]);
    return false;
  }
  if (rank_splits.product() != np) {
    if (rank==0)
      printf("Invalid rank splits. np=%d, product of splits=%d\n",
             np, rank_splits.product());
    return false;
  }
  argno++;

  if (1 != sscanf(argv[argno], "%d", &grid_count)) {
    if (rank==0)
      printf("Invalid grid count: %s\n", argv[argno]);
    return false;
  }
  argno++;

  return true;
}


void fail() {
  MPI_Finalize();
  exit(1);
}


bool setupData(Params &p, const Triple &rank_splits,
               int stripe_count, int stripe_len) {
  // create a cartesian communicator
  int periodic[3] = {0, 0, 0};
  MPI_Cart_create(MPI_COMM_WORLD, 3, rank_splits.coord, periodic, 1,
                  &p.cart_comm);

  // compute my 1-d and 3-d coordinates in the cartesian communicator
  MPI_Comm_rank(p.cart_comm, &p.cart_rank);
  MPI_Cart_coords(p.cart_comm, p.cart_rank, 3, p.cart_coords.coord);

  // determine which subcube this rank is assigned
  for (int i=0; i < 3; i++) {
    partition(p.global_size[i], rank_splits[i], p.cart_coords[i],
              p.local_starts[i], p.local_counts[i]);
  }

  // debug output
  /*
  for (int i=0; i < np; i++) {
    if (i==rank) {
      printf("[%d] cart_rank=%d cart_coords=%s grid=[%d-%d, %d-%d, %d-%d]\n",
             rank, p.cart_rank, p.cart_coords.str(),
             p.local_starts[0], p.local_starts[0] + p.local_counts[0]-1,
             p.local_starts[1], p.local_starts[1] + p.local_counts[1]-1,
             p.local_starts[2], p.local_starts[2] + p.local_counts[2]-1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  */

  p.data.resize(p.local_counts.product());

  /* fill the data */
  for (int z=0; z < p.local_counts[0]; z++) {
    int z_offset = z * p.local_counts[1] * p.local_counts[2];

    for (int y=0; y < p.local_counts[1]; y++) {
      int y_offset = y * p.local_counts[2];

      for (int x=0; x < p.local_counts[2]; x++) {

        p.data[z_offset + y_offset + x] = rank
          + .01 * (z + p.local_starts[0])
          + .0001 * (y + p.local_starts[1])
          + .000001 * (x + p.local_starts[2]);
      }
    }
  }
  
  /* set the file striping parameters */
  MPI_Info_create(&p.info);
  char tmp_str[50];
  sprintf(tmp_str, "%d", stripe_count);
  MPI_Info_set(p.info, "striping_factor", tmp_str);
  sprintf(tmp_str, "%d", stripe_len);
  MPI_Info_set(p.info, "striping_unit", tmp_str);
  
  return true;
}

void printHelp() {
  fprintf(stderr, "\n"
          "  pc2orio2 [opt] <filename> <global_size> <rank_split> <grid_count>\n"
          "    global_size and rank_split are in the form x,y,z\n"
          "  opt:\n"
          "    -stripe_count <n> : number of file stripes\n"
          "    -stripe_len <n> : length of file stripes, in MiB\n"
          "  Compare performance of HDF5 collective IO vs. MPI-IO direct.\n"
          "\n");
  exit(1);
}


void partition(int total_size, int np, int rank, 
               int &offset, int &size) {
  offset = rank * total_size / np;
  int end = (rank + 1) * total_size / np;
  size = end - offset;
}


void timestamp(const char *name) {
  if (PRINT_TIMESTAMPS && rank == 0) {
    printf("%.6f %s\n", MPI_Wtime() - time0, name);
  }
}


int h5_check_internal(int status, const char *filename, int line_no) {
  if (status < 0)
    printf("[%d] %s:%d status %d\n", rank, filename, line_no, status);
  return status;
}


/*
  Write the data using HDF5. This calls MPI_File_write_at_all() internally.
*/
void writeHDF5File(const char *filename, Params &p) {
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
  H5Pset_fapl_mpio(file_access_properties, comm, p.info);

  remove(filename);

  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, 
                      file_create_properties, file_access_properties);
  if (file_id < 0) {
    fprintf(stderr, "Failed to create %s.\n", filename);
    exit(1);
  }

  H5Pclose(file_access_properties);

  timestamp("HDF5 file created");

  MPI_Barrier(MPI_COMM_WORLD);
  timer_open = MPI_Wtime() - timer_open;

  hsize_t global_size[3], local_starts[3], local_counts[3];
  p.global_size.toHsize(global_size);
  p.local_starts.toHsize(local_starts);
  p.local_counts.toHsize(local_counts);

  timer_write = MPI_Wtime();

  for (int grid_no=0; grid_no < p.grid_count; grid_no++) {
    char dataset_name[40];
    sprintf(dataset_name, "/data.%d", grid_no);

    /* create dataspace */
    dataspace_id = H5Screate_simple(3, global_size, NULL);
    h5check(dataspace_id);

    /* create the dataset */
    dataset_id = H5Dcreate2
      (file_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id,
       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    h5check(dataset_id);

    /* define the subset of data I'll write */
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, local_starts,
                                 NULL, local_counts, NULL);
    h5check(status);
  
    /* define the shape of the memory buffer */
    hsize_t mem_offsets[3] = {0, 0, 0};
    memspace_id = H5Screate_simple(3, local_counts, NULL);
    h5check(memspace_id);

    status = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mem_offsets,
                                 NULL, local_counts, NULL);
    h5check(status);

    /* write content */
    hid_t io_prop = H5Pcreate(H5P_DATASET_XFER);
    status = H5Pset_dxpl_mpio(io_prop, H5FD_MPIO_COLLECTIVE);
    h5check(status);

    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,
                      io_prop, /* p.data.data() */ &p.data[0] );
    h5check(status);
    H5Pclose(io_prop);
    h5check(H5Sclose(memspace_id));
    h5check(H5Dclose(dataset_id));
    h5check(H5Sclose(dataspace_id));

    timestamp("dataset written");
  }

  MPI_Barrier(MPI_COMM_WORLD);
  timer_write = MPI_Wtime() - timer_write;

  /* close the file */  
  timer_close = MPI_Wtime();
  h5check(H5Fclose(file_id));
  MPI_Barrier(MPI_COMM_WORLD);
  timer_close = MPI_Wtime() - timer_close;

  if (rank == 0) {
    int64_t bytes = p.grid_count * p.global_size.product() * sizeof(double);
    printf("HDF5 write %s\n"
           "  Open %.3fs, write %.3fs, close %.3fs. Write %.1f MiB/s\n",
           filename,
           timer_open, timer_write, timer_close,
           (bytes / (timer_write * (1<<20))));
  }
}


/*
  Use MPI-IO directly to write the file.
  If use_vector_type use an MPI_Type_vector for the memory array, otherwise
  an MPI_Type_contiguous.
*/
void writeMPIIOFile(const char *filename, Params &p, bool use_vector_type) {
  int rank, status;
  MPI_File file;
  MPI_Datatype file_type, memory_type;
  MPI_Comm comm = MPI_COMM_WORLD;
  double timer_open, timer_write, timer_close;

  MPI_Comm_rank(comm, &rank);

  timer_open = MPI_Wtime();

  MPI_Datatype element_type = MPI_DOUBLE;

  /* define a type for the full array on disk */
  MPI_Type_create_subarray(3, p.global_size.array(), p.local_counts.array(),
                           p.local_starts.array(),
                           MPI_ORDER_C, element_type, &file_type);

  MPI_Type_commit(&file_type);

  /* define a type for the subarray in memory */
  int element_count = p.local_counts.product();
  if (use_vector_type) {
    /* this is what HDF5 uses (H5Smpio.c in H5S_mpio_hyper_type()) */
    MPI_Type_vector(1, element_count, 1, element_type, &memory_type);
  } else {
    /* this is equivalent and more efficient in some cases */
    MPI_Type_contiguous(element_count, element_type, &memory_type);
  }

  MPI_Type_commit(&memory_type);

  /* open the file */
  remove(filename);
  status = MPI_File_open(comm, filename, MPI_MODE_RDWR | MPI_MODE_CREATE,
                         p.info, &file);
  if (status != MPI_SUCCESS) {
    printf("Failed to open \"%s\". Cannot do MPI-IO test.\n", filename);
    return;
  }

  /* set the file view */
  MPI_File_set_view(file, 0, file_type, file_type, "native", p.info);

  MPI_Barrier(comm);
  timer_open = MPI_Wtime() - timer_open;
  timestamp("MPI-IO file created");

  /* write the data */
  MPI_Status status2;
  int count_written;
  timer_write = MPI_Wtime();
  for (int grid_no=0; grid_no < p.grid_count; grid_no++) {
    status = MPI_File_write_at_all(file, grid_no, /* p.data.data() */ &p.data[0], 1,
                                   memory_type, &status2);
    if (status != MPI_SUCCESS) {
      printf("[%d] MPI_File_write failed, error = %d\n", rank, status);
      return;
    }
    MPI_Get_count(&status2, memory_type, &count_written);
    if (count_written != 1)
      printf("[%d] MPI_File_write_at_all count=%d\n", rank, count_written);
    timestamp("grid written");
  }
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
  if (rank == 0) {
    int64_t bytes = p.grid_count * p.global_size.product() * sizeof(double);
    printf("MPI-IO write %s %s\n"
           "  Open %.3fs, write %.3fs, close %.3fs. Write %.1f MiB/s\n",
           use_vector_type ? "vector" : "contiguous", filename,
           timer_open, timer_write, timer_close,
           (bytes / (timer_write * (1<<20))));
  }
}
