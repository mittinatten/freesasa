/*
  Copyright Simon Mitternacht 2013-2015.

  This file is part of FreeSASA.

  FreeSASA is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  FreeSASA is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with FreeSASA.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdarg.h>
#include <stdio.h>
#include "freesasa.h"
#include "util.h"

#ifdef PACKAGE_NAME
const char *freesasa_name = PACKAGE_NAME;
#else
const char *freesasa_name = "freesasa";
#endif

static void
freesasa_err_impl(int err,
                  const char *format,
                  va_list arg)
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
freesasa_mem_fail(const char* func, const char* file, int line)
{
    return freesasa_fail("%s (%s:%d): memory allocation failure.",func,file,line);
}

