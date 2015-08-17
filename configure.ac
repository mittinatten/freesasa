#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([freesasa], [0.2.1.3])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/freesasa.c])
AC_CONFIG_HEADERS([config.h])
m4_include([m4/ax_compare_version.m4])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([pthread], [pthread_create])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h unistd.h pthread.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_C_RESTRICT

# use C99
AC_PROG_CC_C99

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday sqrt strerror])
AC_CHECK_FUNCS([mkdir])

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
AC_CHECK_PROG(DOXYGEN, doxygen, doxygen)

#esyscmd(doxygen --version)]
AM_CONDITIONAL([HAVE_LATEX], [test -n "$LATEX" && test -n "$DVIPS"])
AM_COND_IF([HAVE_LATEX],[],[echo "LaTeX not available, .dvi and .ps documentation cannot be built"])
AM_CONDITIONAL([HAVE_PDF],[test -n "$EPSTOPDF" && test -n "$PDFLATEX"])
AM_COND_IF([HAVE_PDF],[],[.pdf documentation cannot be built])
AM_CONDITIONAL([HAVE_DOXYGEN],[test -n "$DOXYGEN"])
AM_COND_IF([HAVE_DOXYGEN],[AC_CONFIG_FILES([doc/Doxyfile])],[echo "Doxygen not available, html documentation cannot be built"])
if test -n "$DOXYGEN"; then
   DOXYGEN_V=`doxygen --version`
   DOXYGEN_REQ=1.8.8
   AX_COMPARE_VERSION([$DOXYGEN_V],[lt],[$DOXYGEN_REQ],
                      [echo "Warning: Doxygen should preferrably be of version $DOXYGEN_REQ or higher (found $DOXYGEN_V)"])
fi

# Checks for testing
PKG_CHECK_MODULES([CHECK],[check >= 0.9.4],[HAVE_CHECK=1],[HAVE_CHECK=0])
AC_ARG_ENABLE([gcov],
  [AS_HELP_STRING([--enable-gcov],
    [use Gcov to test the test suite (default is no)])],
    [],
    [enable_gcov=no])
AM_CONDITIONAL([COND_GCOV],[test '!' "$enable_gcov" = no])
AM_CONDITIONAL([USE_CHECK],[test "$HAVE_CHECK" -eq 1])

AC_OUTPUT
