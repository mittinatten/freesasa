#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([sasalib], [0.5])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/sasalib.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [sqrt])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [pthread_create])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h unistd.h pthread.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T

# use C99
AC_PROG_CC_C99

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday sqrt strerror])

AC_CONFIG_FILES([Makefile
		 src/Makefile
	         tests/Makefile
                 doc/Makefile
                 doc/fig/Makefile])

# Checks for building manual
AC_CHECK_PROG(LATEX, latex, latex)
AC_CHECK_PROG(DVIPS, dvips, dvips)
AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
AC_CHECK_PROG(EPSTOPDF, epstopdf, epstopdf)

AM_CONDITIONAL([HAVE_LATEX], [test -n "$LATEX" && test -n "$DVIPS"])
AM_COND_IF([HAVE_LATEX],[echo "LaTeX available"],[echo "LaTeX not available, .dvi and .ps documentation cannot be built"])
AM_CONDITIONAL([HAVE_PDF],[test -n "$EPSTOPDF" && test -n "$PDFLATEX"])
AM_COND_IF([HAVE_PDF],[],[.pdf documentation cannot be built])

# Checks for testing
PKG_CHECK_MODULES([CHECK],[check >= 0.9.4])

AC_OUTPUT
