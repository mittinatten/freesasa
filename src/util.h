/*
  Copyright Simon Mitternacht 2013-2016.

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

#ifndef FREESASA_UTIL_H
#define FREESASA_UTIL_H

/**
    @file 
    @author Simon Mitternacht

    Some utility functions and data structures mainly for error
    reporting.
 */

//! Shortcut for memory error generation
#define mem_fail() freesasa_mem_fail(__func__,__FILE__,__LINE__) 

#define fail_msg(msg) freesasa_fail_wloc(__func__,__FILE__,__LINE__,msg)

/**
    Holds interval in a file, to be initalized with ftell() and used
    with fseek().
 */
struct file_interval {
    long begin; //!< Position of beginning of interval
    long end; //!< Position of end of interval
};

//! The name of the library, to be used in error messages and logging
extern const char *freesasa_name;

/**
    Print failure message using format string and arguments.

    @param format Format string
    @return ::FREESASA_FAIL
*/
int 
freesasa_fail(const char *format,...);

/**
    Print warning message using format string and arguments.

    @param format Format string
    @return ::FREESASA_WARN
 */
int
freesasa_warn(const char *format,...);

/**
    Print warning message using function, file and line-number. 

    Usually used via the macro mem_fail().

    @param func Name of the function that failed
    @param file The file where the function is.
    @param line Line number for the error.
    @return ::FREESASA_FAIL
 */
int
freesasa_mem_fail(const char* func,
                  const char* file,
                  int line);

int
freesasa_fail_wloc(const char* func,
                   const char* file,
                   int line,
                   const char *msg);
                

#endif /* FREESASA_UTIL_H */
