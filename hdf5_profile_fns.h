/* 

tau_gen_wrapper minimal.h /usr/local/hdf5/lib/libhdf5.a 

 */

/* #include <stdio.h> */


#define H5_DLL
typedef int herr_t;
typedef char hbool_t;
typedef unsigned long size_t;


H5_DLL herr_t H5open(void);
H5_DLL herr_t H5close(void);
H5_DLL herr_t H5dont_atexit(void);
H5_DLL herr_t H5garbage_collect(void);
H5_DLL herr_t H5set_free_list_limits (int reg_global_lim, int reg_list_lim,
                int arr_global_lim, int arr_list_lim, int blk_global_lim,
                int blk_list_lim);
H5_DLL herr_t H5get_libversion(unsigned *majnum, unsigned *minnum,
				unsigned *relnum);
H5_DLL herr_t H5check_version(unsigned majnum, unsigned minnum,
			       unsigned relnum);
H5_DLL herr_t H5is_library_threadsafe(hbool_t *is_ts);
H5_DLL herr_t H5free_memory(void *mem);
H5_DLL void *H5allocate_memory(size_t size, hbool_t clear);
H5_DLL void *H5resize_memory(void *mem, size_t size);



#include <stdio.h>
int main() {
  size_t x = 0;
  printf("sizeof(size_t)=%d, size_t unsigned: %s \n",
         (int)sizeof(size_t),
         (x-1) > 0 ? "yes" : "no");
  printf("sizeof(bool)=%d\n", (int)sizeof(bool));
  printf("false=%d, true=%d\n", (int)false, (int)true);
  return 0;
}

