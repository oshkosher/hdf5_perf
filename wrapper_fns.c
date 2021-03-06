/* Wrappers of MPI functions to output debugging data. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

extern int rank, np;

#ifdef CRAY_MPI
#define CRAY_CONST const
#else
#define CRAY_CONST
#endif


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


/* Format a string decoding of the "amode" argument to MPI_File_open. */
static char *formatAmode(char buf[174], int amode) {
  /* concatenation of all 9 possible modes with 3-byte separators
     = 173 chars */
  assert(sizeof(amode_flags) / sizeof(amode_flags[0]) == 9);
  char *p = buf;
  *p = 0;
  int i;
  for (i=0; i < 9; i++) {
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


/* Print all the key/value pairs in an MPI_Info object.
   They are output one per line with 'prefix' at the beginning of each line. */
static void printMPIInfo(MPI_Info info, const char *prefix) {
  char key[MPI_MAX_INFO_KEY + 1], *value;
  int value_buf_len = 256;

  if (info == MPI_INFO_NULL) {
    printf("%s  MPI_INFO_NULL\n", prefix);
    return;
  }

  value = (char*) malloc(value_buf_len);
  assert(value);

  int keyno, nkeys;
  MPI_Info_get_nkeys(info, &nkeys);
  for (keyno=0; keyno < nkeys; keyno++) {
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


/* Given an MPI_File, extract the MPI_Info for it and print all
   the key/value pairs with 'prefix' at the beginning of each line
   of output. */
void printFileInfo(MPI_File f, const char *prefix) {
  MPI_Info info;

  MPI_File_get_info(f, &info);

  printMPIInfo(info, prefix);
  
  MPI_Info_free(&info);
}


/* Return a string representing this datatype. If it's a primitive type, return
   the name (MPI_DOUBLE, MPI_INT, ...). Otherwise return a hex string.
   The string is written out 'outstr', and &outstr[0] is returned.
 */
static char* mpiTypeToStr(char outstr[22], MPI_Datatype dt) {
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


/* simple string buffer class for building mpiTypeSignature() result */
typedef struct {
  char *str;
  int len, cap; /* len doesn't include nul terminator */
} StringBuf;

static void sbuf_init(StringBuf *sbuf) {
  sbuf->len = 0;
  sbuf->cap = 1024;
  sbuf->str = (char*) malloc(sbuf->cap);
  sbuf->str[0] = 0;
}

/* make sure the buffer can take an additional 'needed' bytes.
   Return a pointer to the end of the string: str+len */
static char* sbuf_reserve(StringBuf *sbuf, int needed) {
  if (sbuf->len + needed + 1 > sbuf->cap) {
    if (sbuf->cap + needed > sbuf->cap*2)
      sbuf->cap += needed;
    else
      sbuf->cap *= 2;
    sbuf->str = (char*) realloc(sbuf->str, sbuf->cap);
  }
  return sbuf->str + sbuf->len;
}

static char* sbuf_append(StringBuf *sbuf, const char *str) {
  int add_len = strlen(str);
  char *end = sbuf_reserve(sbuf, add_len);
  strcpy(end, str);
  sbuf->len += add_len;
  return sbuf->str + sbuf->len;
}

static char* sbuf_append_int(StringBuf *sbuf, long n) {
  char buf[30];
  sprintf(buf, "%ld", n);
  return sbuf_append(sbuf, buf);
}

/* Returns a string summarizing how the given datatype was constructed.
   The caller must free the result.
 */
char* mpiTypeSignature(MPI_Datatype datatype, int *is_primitive) {
  int i, num_ints, num_addrs, num_dt, combiner;
  MPI_Type_get_envelope(datatype, &num_ints, &num_addrs, &num_dt, &combiner);
  if (combiner == MPI_COMBINER_NAMED) {
    if (is_primitive) *is_primitive = 1;
    char *dt_name = (char*) malloc(30);
    return mpiTypeToStr(dt_name, datatype);
  }

  if (is_primitive) *is_primitive = 0;
  
  StringBuf sbuf;
  sbuf_init(&sbuf);

  const char *combiner_name;
  switch (combiner) {
  case MPI_COMBINER_DUP: combiner_name = "MPI_Type_dup"; break;
  case MPI_COMBINER_CONTIGUOUS: combiner_name = "MPI_Type_contiguous"; break;
  case MPI_COMBINER_VECTOR: combiner_name = "MPI_Type_vector"; break;
  case MPI_COMBINER_HVECTOR: combiner_name = "MPI_Type_create_hvector"; break;
  case MPI_COMBINER_INDEXED: combiner_name = "MPI_Type_indexed"; break;
  case MPI_COMBINER_HINDEXED: combiner_name = "MPI_Type_create_hindexed"; break;
  case MPI_COMBINER_INDEXED_BLOCK: combiner_name = "MPI_Type_create_indexed_block"; break;
  case MPI_COMBINER_HINDEXED_BLOCK: combiner_name = "MPI_Type_create_hindexed_block"; break;
  case MPI_COMBINER_STRUCT: combiner_name = "MPI_Type_create_struct"; break;
  case MPI_COMBINER_SUBARRAY: combiner_name = "MPI_Type_create_subarray"; break;
  case MPI_COMBINER_DARRAY: combiner_name = "MPI_Type_create_darray"; break;
  case MPI_COMBINER_F90_REAL: combiner_name = "MPI_Type_create_F90_real"; break;
  case MPI_COMBINER_F90_COMPLEX: combiner_name = "MPI_Type_create_F90_complex"; break;
  case MPI_COMBINER_F90_INTEGER: combiner_name = "MPI_Type_create_F90_integer"; break;
  case MPI_COMBINER_RESIZED: combiner_name = "MPI_Type_create_resized"; break;
  default: combiner_name = "<unknown>"; break;
  }

  int *ints = (int*) malloc(sizeof(int) * num_ints);
  MPI_Aint *addrs = (MPI_Aint*) malloc(sizeof(MPI_Aint) * num_addrs);
  MPI_Datatype *datatypes = (MPI_Datatype*) malloc(sizeof(MPI_Datatype) * num_dt);
  MPI_Type_get_contents(datatype, num_ints, num_addrs, num_dt, ints,
                        addrs, datatypes);
  
  sbuf_append(&sbuf, combiner_name);
  sbuf_append(&sbuf, "(ints={");
  for (i=0; i < num_ints; i++) {
    if (i > 0) sbuf_append(&sbuf, ", ");
    sbuf_append_int(&sbuf, ints[i]);
  }
  sbuf_append(&sbuf, "}, addrs={");
  for (i=0; i < num_addrs; i++) {
    if (i > 0) sbuf_append(&sbuf, ", ");
    sbuf_append_int(&sbuf, addrs[i]);
  }
  sbuf_append(&sbuf, "}, types={");
  for (i=0; i < num_dt; i++) {
    if (i > 0) sbuf_append(&sbuf, ", ");
    int is_prim;
    char *dt_str = mpiTypeSignature(datatypes[i], &is_prim);
    sbuf_append(&sbuf, dt_str);
    free(dt_str);
    if (!is_prim) MPI_Type_free(&datatypes[i]);
  }
  sbuf_append(&sbuf, "})");
  return sbuf.str;
}


int MPI_File_open
  (MPI_Comm comm, const char *filename, int amode, MPI_Info info,
   MPI_File *fh) {
  const char *prefix = "[0]  ";

  /* amode |= MPI_MODE_UNIQUE_OPEN; */

  /*
  if (rank==0) {
    printf("[%d] MPI_File_open %s info:\n", rank, filename);
    printMPIInfo(info, prefix);
  }
  */

  int result = PMPI_File_open(comm, filename, amode, info, fh);
  
  if (rank==0) {
    char amode_str[200];
    printf("[%d] MPI_File_open(\"%s\", amode=%s, fh=%p)\n",
           rank, filename, formatAmode(amode_str, amode), (void*)*fh);
    /* printFileInfo(*fh, prefix); */
  }
  
  return result;
}


int MPI_File_sync(MPI_File fh) {
  printf("[%d] MPI_file_sync(%p)\n", rank, (void*)fh);
  return PMPI_File_sync(fh);
}


int MPI_file_seek(MPI_File fh, MPI_Offset offset, int whence) {
  printf("[%d] MPI_File_seek(%p, %ld, %s)\n",
         rank, (void*)fh, (long)offset,
         whence == MPI_SEEK_SET ? "SET" :
         whence == MPI_SEEK_CUR ? "CUR" : 
         whence == MPI_SEEK_END ? "END" : "<unknown>");
  return MPI_File_seek(fh, offset, whence);
}


int MPI_File_set_view
  (MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,
   const char *datarep, MPI_Info info) {
  const char *prefix = "[0]  ";
  
  if (rank==0) {
    // char buf1[50], buf2[50];
    char *buf1 = mpiTypeSignature(etype, NULL), *buf2 = mpiTypeSignature(filetype, NULL);
    printf("[%d] MPI_File_set_view fh=%p disp=%ld, etype=%s, filetype=%s, datarep=%s\n",
           // rank, (void*)fh, (long)disp, mpiTypeToStr(buf1, etype), mpiTypeToStr(buf2, filetype), datarep);
           rank, (void*)fh, (long)disp, buf1, buf2, datarep);
    /* printMPIInfo(info, prefix); */
    free(buf1);
    free(buf2);
  }

  return PMPI_File_set_view(fh, disp, etype, filetype, datarep, info);
}


int MPI_File_set_size(MPI_File fh, MPI_Offset size) {
  if (rank==0)
    printf("[%d] MPI_File_set_size(%p, %ld)\n", rank, (void*)fh, (long)size);
  
  return PMPI_File_set_size(fh, size);
}


int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void * buf, int count, MPI_Datatype datatype, MPI_Status * status) {
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_File_write_at(fh=%p, offset=%ld, count=%d, type=%s)\n",
           (void*)fh, (long)offset, count, mpiTypeToStr(buf, datatype));
  }
  return PMPI_File_write_at(fh, offset, buf, count, datatype, status);
}

int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void * buf, int count, MPI_Datatype datatype, MPI_Status * status) {
  if (rank==0) {
    char *buf = mpiTypeSignature(datatype, NULL);
    // printf("[0] MPI_File_write_at_all(file=%p, offset=%ld, count=%d, type=0x%08x)\n", (void*)fh, (long)offset, count, (int)datatype);
    printf("[0] MPI_File_write_at_all(file=%p, offset=%ld, count=%d, type=%s)\n", (void*)fh, (long)offset, count, buf);
  }
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
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_contiguous(count=%d, oldtype=%s, newtype=0x%08x)\n",
           count, mpiTypeToStr(buf, oldtype), *(int*)newtype);
  }
  */
  return result;
}

int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype * newtype) {
  int result = PMPI_Type_vector(count, blocklength, stride, oldtype, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_vector(count=%d, blocklength=%d, stride=%d, oldtype=%s, newtype=0x%08x)\n",
           count, blocklength, stride, mpiTypeToStr(buf, oldtype), *(int*)newtype);
  }
  */
  return result;
}

int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype * newtype) {
  int result = PMPI_Type_hvector(count, blocklength, stride, oldtype, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_hvector(%d, %d, %ld, %s, 0x%08x)\n",
           count, blocklength, stride, mpiTypeToStr(buf, oldtype), *(int*)newtype);
  }
  */
  return result;
}

int MPI_Type_indexed(int count, const int * array_of_blocklengths, const int * array_of_displacements, MPI_Datatype oldtype, MPI_Datatype * newtype) {
  int result = PMPI_Type_indexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_indexed(%d, ..., ..., %s, 0x%08x)\n",
           count, mpiTypeToStr(buf, oldtype), *(int*)newtype);
  }
  */
  return result;
}

int MPI_Type_hindexed(int count, CRAY_CONST int * array_of_blocklengths, CRAY_CONST MPI_Aint * array_of_displacements, MPI_Datatype oldtype, MPI_Datatype * newtype) {
  int result = PMPI_Type_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_hindexed(%d, ..., ..., %s, 0x%08x)\n",
           count, mpiTypeToStr(buf, oldtype), *(int*)newtype);
  }
  */
  return result;
}

int MPI_Type_struct(int count, CRAY_CONST int * array_of_blocklengths, CRAY_CONST MPI_Aint * array_of_displacements, CRAY_CONST MPI_Datatype * array_of_types, MPI_Datatype * newtype) {
  int result = PMPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_struct(%d, ..., ..., 0x%08x)\n",
           count, *(int*)newtype);
  }
  */
  return result;
}


int MPI_Type_create_darray(int size , int rankarg , int ndims , const int array_of_gsizes [], const int array_of_distribs [], const int array_of_dargs [], const int array_of_psizes [], int order , MPI_Datatype oldtype , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_darray(size, rankarg, ndims, array_of_gsizes, array_of_distribs, array_of_dargs, array_of_psizes, order, oldtype, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_darray(... %s, 0x%08x)\n",
           mpiTypeToStr(buf, oldtype), *(int*)newtype);
  }
  */
  return result;
}

int MPI_Type_create_hindexed(int count , const int array_of_blocklengths [], const MPI_Aint array_of_displacements [], MPI_Datatype oldtype , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_hindexed(count=%d, oldtype=%s, newtype=0x%08x)\n",
           count, mpiTypeToStr(buf, oldtype), *(int*)newtype);
    
    for (int i=0; i < count; i++) {
      if (i > 0) printf(", ");
      printf("%d@%ld", array_of_blocklengths[i], (long)array_of_displacements[i]);
    }
    putchar('\n');

  }
  */
  return result;
}

int MPI_Type_create_hvector(int count , int blocklength , MPI_Aint stride , MPI_Datatype oldtype , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_hvector(count, blocklength, stride, oldtype, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_hvector(... %s, 0x%08x)\n",
           mpiTypeToStr(buf, oldtype), *(int*)newtype);
  }
  */
  return result;
}

int MPI_Type_create_indexed_block(int count , int blocklength , const int array_of_displacements [], MPI_Datatype oldtype , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_indexed_block(count, blocklength, array_of_displacements, oldtype, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_indexed_block(... %s, 0x%08x)\n",
           mpiTypeToStr(buf, oldtype), *(int*)newtype);
  }
  *?
  return result;
}

int MPI_Type_create_hindexed_block(int count , int blocklength , const MPI_Aint array_of_displacements [], MPI_Datatype oldtype , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_hindexed_block(count, blocklength, array_of_displacements, oldtype, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_hindexed_block(... %s, 0x%08x)\n",
           mpiTypeToStr(buf, oldtype), *(int*)newtype);
  }
  */
  return result;
}

int MPI_Type_create_resized(MPI_Datatype oldtype , MPI_Aint lb , MPI_Aint extent , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_resized(oldtype, lb, extent, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_resized(lb=%ld, extent=%ld, oldtype=%s, newtype=0x%08x)\n",
           (long)lb, (long)extent, mpiTypeToStr(buf, oldtype), *(int*)newtype);
  }
  */
  return result;
}

int MPI_Type_create_struct(int count , const int array_of_blocklengths [], const MPI_Aint array_of_displacements [], const MPI_Datatype array_of_types [], MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_struct(count=%d, newtype=0x%08x) {", count, *(int*)newtype);
    for (int i=0; i < count; i++) {
      char buf[50];
      if (i > 0) printf(", ");
      printf("%d %s @ %ld", array_of_blocklengths[i], mpiTypeToStr(buf, array_of_types[i]), (long)array_of_displacements[i]);
    }
    printf("}\n");
  }
  */
  return result;
}

int MPI_Type_create_subarray(int ndims , const int array_of_sizes [], const int array_of_subsizes [], const int array_of_starts [], int order , MPI_Datatype oldtype , MPI_Datatype * newtype ) {
  int result = PMPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, order, oldtype, newtype);
  /*
  if (rank==0) {
    char buf[50];
    printf("[0] MPI_Type_create_subarray(... %s, 0x%08x)\n",
           mpiTypeToStr(buf, oldtype), *(int*)newtype);
  }
  */
  return result;
}

int MPI_Type_commit(MPI_Datatype *datatype) {
  int result = PMPI_Type_commit(datatype);
  /*
  if (rank==0) {
    printf("[0] MPI_Type_commit(0x%08x)\n", *datatype);
  }
  */
  return result;
}


int MPI_Type_free(MPI_Datatype *datatype) {
  /*
  if (rank==0) {
    printf("[0] MPI_Type_free(0x%08x)\n", *datatype);
  }
  */
  return PMPI_Type_free(datatype);
}

