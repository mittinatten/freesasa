#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include "freesasa.h"
#include "freesasa_internal.h"

#ifdef PACKAGE_NAME
const char *freesasa_name = PACKAGE_NAME;
#else
const char *freesasa_name = "freesasa";
#endif

static FILE *errlog = NULL;

struct file_range
freesasa_whole_file(FILE* file)
{
    assert(file);
    struct file_range range;
    rewind(file);
    range.begin = ftell(file);
    fseek(file,0,SEEK_END);
    range.end = ftell(file);
    rewind(file);
    assert(range.begin <= range.end);
    return range;
}

static void
freesasa_err_impl(int err,
                  const char *format,
                  va_list arg)
{
    FILE *fp = stderr;
    if (errlog != NULL) fp = errlog;

    fprintf(fp, "%s: ", freesasa_name);
    switch (err) {
    case FREESASA_FAIL: fputs("error: ", fp); break;
    case FREESASA_WARN: fputs("warning: ", fp); break;
    default: break;
    }
    vfprintf(fp, format, arg);
    va_end(arg);
    fputc('\n', fp);
    fflush(fp);
}

int
freesasa_fail(const char *format,...)
{
    va_list arg;
    if (freesasa_get_verbosity() == FREESASA_V_SILENT) return FREESASA_FAIL;
    va_start(arg, format);
    freesasa_err_impl(FREESASA_FAIL,format,arg);
    va_end(arg);
    return FREESASA_FAIL;
}

int freesasa_warn(const char *format,...)
{
    va_list arg;
    int v = freesasa_get_verbosity();
    if (v == FREESASA_V_NOWARNINGS || v == FREESASA_V_SILENT) return FREESASA_WARN;
    va_start(arg, format);
    freesasa_err_impl(FREESASA_WARN,format,arg);
    va_end(arg);
    return FREESASA_WARN;
}

int
freesasa_fail_wloc(const char* func,
                   const char* file,
                   int line,
                   const char *msg) 
{
    return freesasa_fail("in %s() (%s:%d): %s",func,file,line,msg);
}

int
freesasa_mem_fail(const char* func, const char* file, int line)
{
    return freesasa_fail_wloc(func,file,line,"Out of memory.");
}

const char*
freesasa_thread_error(int error_code)
{
    switch (error_code) {
    case EDEADLK: return "deadlock (EDEADLK)";
    case EINVAL: return "(EINVAL)";
    case ESRCH: return "no matching thread (ESRCH)";
    case EAGAIN: return "no resources to create thread (EAGAIN)";
    }
    return "Unknown thread error";
}

void
freesasa_set_err_out(FILE *fp)
{
    assert(fp);
    errlog = fp;
}

FILE *
freesasa_get_err_out()
{
    return errlog;
}
