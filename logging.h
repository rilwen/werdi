#ifndef _LOGGING_H
#define _LOGGING_H

#include <stdio.h>

/* Return the pointer to the log file */
FILE* get_log_file(void);

/* Start logging */
void start_loggin(void);

/* Stop logging */
void stop_logging(void);

#endif
