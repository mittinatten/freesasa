#ifndef UTIL_H
#define UTIL_H

#define mem_fail() freesasa_mem_fail(__func__,__FILE__,__LINE__) 

struct file_interval {long begin; long end;};

extern const char *freesasa_name;

int freesasa_fail(const char *format,...);

int freesasa_warn(const char *format,...);

int freesasa_mem_fail(const char* func, const char* file, int line);

#endif /* UTIL_H */
