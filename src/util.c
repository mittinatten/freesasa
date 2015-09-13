#include <stdarg.h>
#include <stdio.h>
#include "freesasa.h"
#include "util.h"

extern const char* freesasa_name;

static void freesasa_err_impl(int err, const char *format, va_list arg)
{
    fprintf(stderr, "%s: ", freesasa_name);
    switch (err) {
    case FREESASA_FAIL: fputs("error: ", stderr); break;
    case FREESASA_WARN: fputs("warning: ", stderr); break;
    default: break;
    }
    vfprintf(stderr, format, arg);
    va_end(arg);
    fputc('\n', stderr);
    fflush(stderr);
}

int freesasa_fail(const char *format,...)
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

int freesasa_mem_fail(const char* func, const char* file, int line)
{
    return freesasa_fail("%s (%s:%d): memory allocation failure.",func,file,line);
}

