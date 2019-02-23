#ifndef __IO_WRAPPERS_H__
#define __IO_WRAPPERS_H__

#include <stdio.h>


#ifdef __cplusplus
extern "C" {
#endif

void io_wrappers_reset();
void io_wrappers_report(FILE *outf, double *io_time, double *comm_time);

#ifdef __cplusplus
}
#endif

#endif /* __IO_WRAPPERS_H__ */
