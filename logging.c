#include <stdio.h>

FILE* global_log_file = NULL;

FILE* get_log_file(void)
{
	return global_log_file;
}

void start_loggin(void)
{
	if (global_log_file == NULL)
	{
		global_log_file = fopen("werdi_log.txt", "a");
	}
}

void stop_logging(void)
{
	if (global_log_file != NULL)
	{
		fclose(global_log_file);
	}
}
