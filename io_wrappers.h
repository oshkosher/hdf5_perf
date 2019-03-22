#ifndef __IO_WRAPPERS_H__
#define __IO_WRAPPERS_H__

#include <stdio.h>

/*
  If the environment variable "IO_WRAPPERS_LOG_PREFIX" is set, it
  is used as the name of the MPE logfile, with ".clog2" appended.
  Otherwise the name "io_wrappers.clog2" is used.
*/

#ifdef __cplusplus
extern "C" {
#endif

void io_wrappers_reset();
void io_wrappers_report(FILE *outf);


#ifdef __cplusplus
}
#endif

#endif /* __IO_WRAPPERS_H__ */
