#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([sasalib], [0.5])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/sasalib.c])
#AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [sqrt])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [pthread_create])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday sqrt strerror])

AC_CONFIG_FILES([Makefile
		 src/Makefile
	         tests/Makefile
                 doc/Makefile
                 doc/fig/Makefile])

AC_CHECK_PROG(LATEX, latex, latex)
AC_CHECK_PROG(DVIPS, dvips, dvips)
AC_CHECK_PROG(PS2PDF, ps2pdf, ps2pdf)

AM_CONDITIONAL([HAVE_LATEX], test -n "$PDFLATEX" && test -n "$DVIPS" && test -n "$PS2PDF")

AC_OUTPUT
